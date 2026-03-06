import numpy as np
import pandas as pd
from PIL import Image as im
import csv
import circle_fit
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from scipy.interpolate import interp1d
import os
import re

import pretty_errors

def rename_files(folder):

    files = os.listdir(folder)
    print(files)

    for file in files:
        print(file)
        if file.endswith('.csv'):
            RunList = file
            RunDict = csv.DictReader(open(folder + '/' + RunList, encoding = 'utf-8-sig'))

    print(RunDict)
    #
    filenum = re.compile(r"\d\d\d\d\d\d\d")
    batch_info_dicts = []

    for dict in RunDict:
        batch_info_dicts.append(dict)

    for file in files:
        if file.endswith('.tiff'):
            match = filenum.search(file)
            print(file)
            ID = match.group()
            for dict in batch_info_dicts:
                print(dict)
                if ID in dict['Filename']:
                    newName = dict['Sample Description'][12:] + '.tiff'
                    print(newName)
                    os.rename(folder + '/' + file, folder + '/' + newName)
                    print(dict['Filename'] + ' renamed to ' + newName)

                else:
                    pass
    return

def process_run_list(path):
    run_dict = csv.DictReader(open(path))
    return run_dict

def get_norm_factor(sample_name, run_dict, distance_mm, pixel_size_mm):
    for sample_dict in run_dict:

        if sample_name == sample_dict['Sample Description'][12:]:
            scalar = 1 / (pixel_size_mm / distance_mm) ** 2 / float(sample_dict['Measurement Time']) / float(
                sample_dict['Transmission']) / float(sample_dict['Io'])
            break

    return scalar

def tiff_to_pixel_grid(tiff_array):
    a = tiff_array
    # Get image shape
    height, width = a.shape

    # Create X, Y grid
    x, y = np.meshgrid(np.arange(width), np.arange(height))
    z = a

    # Flatten to XYZ
    xyz = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

    return xyz

def grid_to_polar(xyz, center):
    x = xyz[:, 0] - center[0]
    y = xyz[:, 1] - center[1]
    z = xyz[:, 2]

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    z = z
    rtz = np.column_stack((r,theta,z))
    return rtz

def polar_to_grid(rtz,center):
    r = rtz[:, 0]
    t= rtz[:, 1]
    z = rtz[:,2]

    x = r * np.cos(t) + center[0]
    y = r * np.sin(t)+center[1]
    z = z

    xyz = np.column_stack((x,y,z))
    return xyz

def match_xyz(xyz, indx):
    loc = np.where((xyz[:, 0] == indx[0]) & (xyz[:, 1] == indx[1]))
    point = xyz[loc][0]
    return point

def mask_theta(agb_xyz, center):

    x0 = center[0]
    y0=center[0]
    mask0 = agb_xyz[:, 2] < 0

    points = agb_xyz[agb_xyz[:, 2] >= 0]


    xi = points[:, 0]-x0
    yi= points[:, 1]-y0
    zi = points[:, 2]

    rtz = grid_to_polar(points,center)
    beamstop = 8
    x,y,y_error, bin_pixel_count = bin_points(rtz[rtz[:,0]>beamstop],1,num_bins = 360)
    theta_fig, theta_ax = plt.subplots()
    theta_ax.set_xlabel('Angle')
    theta_ax.set_ylabel('I (arb)')
    theta_ax.plot(x,y)
    def select_mask(xmin, xmax):
        plt.close(theta_fig)

    span = SpanSelector(theta_ax, select_mask,'horizontal')
    theta_min, theta_max = span.extents
    plt.show()
    # diff_y = np.diff(y)
    # diff_x = np.diff(x)
    # dydx = diff_y/diff_x
    # mean_diff = np.mean(dydx)
    # xmid = 0.5*(x[:-1]+x[1:])
    # bad_diff = []
    # chunk = 1
    # for i in range(0,len(dydx),chunk):
    #     a = dydx[i-chunk:i+chunk]
    #     chunk_avg = np.mean(chunk)
    #     if all(abs(a)<10000) and i>5 and i<len(dydx-chunk):
    #         bad_diff.append(i-chunk)
    #         bad_diff.append(i+chunk)

    # bad_x = xmid[bad_diff[1:]]

    masked_theta_range = span.extents

    return masked_theta_range

def agb_cal(tiff_array):
    array = tiff_array
    xyz = tiff_to_pixel_grid(array)
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]

    def select_ring(array):
        agb_fig = plt.figure('Initial', figsize=(8, 8))
        ax_agb = agb_fig.add_subplot(111)
        ax_agb.imshow(np.log10(array), cmap='terrain')
        ax_agb.set_facecolor('black')
        ax_agb.invert_yaxis()
        ring_point = plt.ginput(1)
        plt.close()

        ring_index = np.int_(np.round(ring_point[0])) - 1
        matched = match_xyz(xyz, ring_index)
        return matched
    def find_ring(start):
        sys.setrecursionlimit(30000)
        visited = set()
        valid = set()
        valid.add(start)

        def neighbors(point):
            offset = ((-1,0),(1,0),(0,-1),(0,1))
            for dx, dy in offset:
                x, y = point[0] + dx, point[1] + dy
                new = match_xyz(xyz, (x,y))
                yield tuple(new)

        def branch(point):
            thresh = point[2] * 0.1
            for neighbor in neighbors(point):
                if neighbor not in visited:
                    visited.add(tuple(neighbor))

                    if neighbor[2]>=thresh:
                        valid.add(tuple(neighbor))
                        branch(tuple(neighbor))
                    else:
                        return
        branch(start)
        return valid


    xyz0 = tuple(select_ring(array))

    fig_ring, ax_ring = plt.subplots()
    fig_ring.suptitle('AgB first ring calibration')
    ax_ring.imshow(np.log10(array))
    ax_ring.invert_yaxis()
    ax_ring.set_xlabel('Horizontal pixels')
    ax_ring.set_ylabel('Vertical pixels')

    ring_points = find_ring(xyz0)

    x0,y0, r_pixel, sig =  circle_fit.hyperLSQ(list(ring_points))

    rp_array = np.array((list(ring_points)))


    ax_ring.plot(rp_array[:,0], rp_array[:,1], marker = '.', markersize = 1, markerfacecolor = 'red', markeredgecolor = 'red', ls = '')
    ax_ring.plot(x0,y0,ls='', marker = '*', markersize =1, markerfacecolor = 'red', markeredgecolor='red')
    plt.show()
    r_mm = 0.172 * r_pixel
    d_AgB = 58.380 * 1E-7
    # angstroms to mm
    distance = r_mm * d_AgB / (1.54 * 1E-7)

    return [x0,y0], distance,r_pixel, rp_array, fig_ring

def bin_points(points,indx,binning = 'linear',num_bins=400,minmax = None):
    bin_var = points[:,indx]
    if minmax is None:
        minmax = [min(bin_var),max(bin_var)]

    if binning == 'linear':
        bin_edges = np.linspace(minmax[0], minmax[1],num_bins)

    elif binning == 'log':
        bin_edges = np.logspace(np.log10(minmax[0]),np.log10(minmax[1]),num_bins)

    bindex = np.digitize(bin_var, bin_edges)

    bin_sums = [0]*len(bin_edges)
    sqr_res_sum = [0]*len(bin_edges)
    bin_pixel_count = [0]*len(bin_edges)

    for i, val in enumerate(bin_var):
        bin = bindex[i]-1
        bin_sums[bin]+=points[i,2]
        bin_pixel_count[bin]+=1

    bin_avgs = np.array(bin_sums)/np.array(bin_pixel_count)

    for i, val in enumerate(bin_var):
        bin = bindex[i]-1
        sqr_res_sum[bin]+= (points[i,2]-bin_avgs[bin])**2

    bin_errors = np.sqrt(np.array(sqr_res_sum)/np.array(bin_pixel_count)**2)

    ## Poisson Error approach
    # bin_errors = np.sqrt(np.array(bin_sums)/np.array(bin_pixel_count)**2)

    return bin_edges, bin_avgs, bin_errors, bin_pixel_count

def interp_y(x0, y0, x):
    interp_func = interp1d(x0, y0, fill_value='extrapolate')

    return interp_func(x)

def bkg_sub_func(data,bkg):
    q = data[0]
    I = data[1]
    e = data[2]

    I_bkg = interp_y(bkg[0],bkg[1],q)
    e_bkg = interp_y(bkg[0],bkg[2], q)




    new_I = I - I_bkg
    new_e = np.sqrt(e**2 - e_bkg**2)
    return np.array([q,new_I,new_e])

def plot_data(data):

    sp_dict = {'xscale':'log','yscale':'log','xlabel':'q ($A^{-1}$)','ylabel':'I (arb)'}
    print(data.keys())
    if 'raw' in data.keys():
        fig_r, ax_r = plt.subplots(subplot_kw=sp_dict)
        fig_r.suptitle('Raw data')

        for sample, a in data['raw'].items():

            ax_r.errorbar(a[0],a[1],a[2], ls = '', marker = 'o', markersize='5', markerfacecolor = 'None', label=sample)
        ax_r.legend()
    if 'scale' in data.keys():
        fig_s, ax_s = plt.subplots(subplot_kw=sp_dict)
        fig_s.suptitle('Normalized data')
        if 'bkg_data' in data.keys():
            ax_s.errorbar(data['bkg_data'][0], data['bkg_data'][1], data['bkg_data'][2], ls = '', marker = 'o',
                        markersize = '5', label = 'background')
        for sample, a in data['scale'].items():
            ax_s.errorbar(a[0],a[1],a[2],ls = '', marker = 'o', markerfacecolor = 'None', markersize = '5', label=sample)
        ax_s.legend()
    if 'bkg_sub' in data.keys():
        fig_b, ax_b = plt.subplots(subplot_kw=sp_dict)
        fig_b.suptitle('Background Subtracted')

        for sample, a in data['bkg_sub'].items():
            ax_b.errorbar(a[0],a[1],a[2], ls='', marker = 'o', markersize = '5', markerfacecolor = 'None', label=sample)
        ax_b.legend()
    plt.show()

def save_data(data, save_folder, figs, exp_params):
    print(save_folder)
    # with open(save_folder, 'data_parameters.txt','r') as f:
    #     for key, value in exp_params.items():
    #         f.write(f'{key}: {value}\n')
    #     f.close()

    agb_figure_path = os.path.join(save_folder, 'agb_figure')
    os.makedirs(agb_figure_path, exist_ok=True)


    figs['AgB_ring'].savefig(os.path.join(agb_figure_path, 'AgB_ring_calibration.png'))
    figs['mask'].savefig(os.path.join(agb_figure_path, 'mask.png'))

    csv_path = os.path.join(save_folder, 'csv')
    os.makedirs(csv_path, exist_ok=True)
    
    npy_path = os.path.join(save_folder, 'npy')
    os.makedirs(npy_path, exist_ok=True)

    def treat_csv(name, a, parent):
        df = pd.DataFrame(a.transpose(), columns=['q (1/A)', 'Intensity (arb)', 'Error (arb)'])
        save_path = os.path.join(parent, name + '.csv')
        df.to_csv(save_path, index=False)
        return

    for type in ['raw', 'scale', 'bkg_sub']:
        if type in data.keys():
            parent_csv = os.path.join(csv_path,type)
            os.makedirs(parent_csv, exist_ok=True)

            parent_npy = os.path.join(npy_path,type)
            os.makedirs(parent_npy, exist_ok=True)
            for sample, a in data[type].items():
                treat_csv(sample+'_'+type,a,parent_csv)
                npy_name = os.path.join(parent_npy,sample)
                np.save(npy_name, a)


    

    # ## csv
    # for parent in ['csv','npy']:
    #     parent_path = os.path.join(save_folder, parent)
    #     os.makedirs(parent, exist_ok=True)
    #     raw_parent = os.path.join(parent_path, 'raw')
    #     os.makedirs(raw_parent, exist_ok=True)
    #     
    #     for sample, a in data['raw'].items():
    #         if parent =='csv':
    #             treat_csv(sample,a, raw_parent)
    #         if parent == 'npy'
    #             np.save(os.path.join(raw_parent, sample), a)
    #     
    #     if 'scale' in data.keys():
    #         scale_parent = os.path.join(parent_path, 'normalized')
    #         
    #         for sample, a in data['scale'].items():
    #             if parent == 'csv':
    #                 treat_csv(sample,a,scale_parent)
    #             if parent == 'npy':
    #                 np.save(os.path.join(scale_parent, sample), a)
    #     
    #     if 'bkg_sub' in data.keys():
    #         bkg_sub_parent = os.path.join(parent_path,'bkg_sub')
    #         
    #         for sample, a in data['bkg_sub'].items():
    #             if parent =='csv':
    #                 treat_csv(sample,a,bkg_sub_parent)
    #             if parent == 'npy':
    #                 np.save(os.path.join(bkg_sub_parent,sample),a)

def run_folder(folder,scale=True,bkg = None, plot=True, save = True, center = None,apply_mask = True,mask = None,
               binning = {'var':'q','style':'log','num_bins':400},pixel_size = 0.172, rename = False):
    if rename:
        rename_files(folder)

    files = os.listdir(folder)
    files.sort()
    tiff_a = {}
    cal_info = {}
    figs = {}
    exp_params = {}
    for f in files:
        path = os.path.join(folder,f)


        if f.endswith('.tiff'):
            name = os.path.basename(f).removesuffix('.tiff')
            if 'AgB' in f or 'agb' in f or 'Agb' in f and '.tiff' in f:
                agb_a = np.array(im.open(path))
            elif f.endswith('tiff'):
                tiff_a[name] = np.array(im.open(path))
        elif f.endswith('.csv'):
            run_list_path = os.path.join(folder,f)
            exp_params['run_list'] = run_list_path

    if center == None:
        center, distance, r_pixel, ring_points, figs['AgB_ring'] = agb_cal(agb_a)
        cal_info['center'] = center
        cal_info['distance'] = distance
        exp_params['pixel_center'] = center
        exp_params['detector_distance (mm)']=distance
    # center = [267.46,359.05]
    # distance = 1079.64
    # r_pixel = 165.58

    data = {}

    raw = {}
    scaled = {}
    bkg_sub = {}
    if apply_mask == True:
        agb_xyz = tiff_to_pixel_grid(agb_a)
        masked_theta_range = mask_theta(agb_xyz,center)
        dead_pix = agb_xyz[agb_xyz[:, 2] < 0]
        agb_xyz = agb_xyz[agb_xyz[:, 2] >= 0]
        agb_rtz = grid_to_polar(agb_xyz, center)

        bp = agb_rtz[agb_rtz[:, 0] < 8]
        agb_rtz = agb_rtz[agb_rtz[:, 0] > 8]

        bp = polar_to_grid(bp, center)

        theta_mask = agb_rtz[(agb_rtz[:, 1] > masked_theta_range[0]) & (agb_rtz[:, 1] < masked_theta_range[1])]
        theta_mask = polar_to_grid(theta_mask,center)


        fig_m, ax_m = plt.subplots()
        ax_m.imshow(np.log10(agb_a), cmap='terrain')
        ax_m.plot(bp[:,0],bp[:,1], ls = '', marker = 's', markersize =1, color = 'black')
        ax_m.plot(dead_pix[:,0],dead_pix[:,1], ls = '', marker = 's', markersize =1, color = 'black')
        ax_m.plot(theta_mask[:,0],theta_mask[:,1], ls = '', marker = 's', markersize =1, color = 'black')
        figs['mask'] = fig_m

    scalars = []
    for sample, a in tiff_a.items():
        xyz = tiff_to_pixel_grid(a)

        if apply_mask:
            xyz = xyz[xyz[:,2]>=0]
            rtz = grid_to_polar(xyz,center)

            rtz = rtz[rtz[:,0]>8]

            rtz = rtz[(rtz[:,1]>masked_theta_range[0]) | (rtz[:,1]<masked_theta_range[1])]


        else:
            rtz = grid_to_polar(xyz,center)


        if binning['var'] == 'q':
            two_theta = np.arctan(rtz[:, 0] * 0.172 / distance)
            q = 4*np.pi*np.sin(two_theta/2)/1.54
            qtz = np.zeros_like(rtz)
            qtz[:,0] = q
            qtz[:,[1,2]] = rtz[:,[1,2]]
            q_edges, I_bins, I_bin_errors, num_points = bin_points(qtz, 0, binning = binning['style'],
                                                                   num_bins= binning['num_bins'])
        elif binning['var'] == 'r':
            r_edges, I_bins, I_bin_errors, num_points = bin_points(rtz,0,binning = binning['style'],
                                                                 num_bins = binning['num_bins'])
            two_theta = np.arctan(r_edges * 0.172 / distance)
            q_edges = 4*np.pi*np.sin(two_theta/2)/1.54

        raw[sample] = np.array([q_edges,I_bins,I_bin_errors])
        df = pd.DataFrame(raw[sample].transpose())

        if scale:

            run_dict = process_run_list(run_list_path)
            norm_factor = get_norm_factor(sample, run_dict, distance, pixel_size)
            scaled[sample]=np.array([q_edges, raw[sample][1] * norm_factor, raw[sample][2] * norm_factor])
            scalars.append([sample, norm_factor])
        if bkg is not None:
            cal_info['background'] = bkg_path
            bkg_data = np.load(bkg)
            bkg_sub[sample]=bkg_sub_func(scaled[sample],bkg_data)



    data['raw'] = raw
    if scale:
        data['scale'] = scaled
        exp_params['scalar'] = scalars
    if bkg is not None:
        data['bkg_data'] = bkg_data
        data['bkg_sub'] = bkg_sub
        exp_params['background_file']= bkg_path
    if save:
        save_data(data,folder,figs, exp_params)

    if plot:
        plot_data(data)

if __name__ == '__main__':
    folder = '/Users/willbrackett/Downloads/SAXS Data'
    bkg_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/Backgrounds/dmf_bkg_60mins_norm.npy'
    run_folder(folder, bkg = bkg_path,rename = False)

