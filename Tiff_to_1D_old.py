import pandas as pd
from PIL import Image as im
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sys
import os
import circle_fit
import re
import csv
import datetime
from scipy.interpolate import interp1d

import tkinter as tk
from tkinter import filedialog

def rename_files(folder):

    files = os.listdir(folder)

    for file in files:
        print(file)
        if file.endswith('.csv'):
            RunList = file
            RunDict = csv.DictReader(open(folder + '/' + RunList))

    print(RunDict)
    #
    filenum = re.compile(r"\d\d\d\d\d\d\d")
    batch_info_dicts = []

    for dict in RunDict:
        batch_info_dicts.append(dict)

    for file in files:
        if file.endswith('.tiff'):
            match = filenum.search(file)
            ID = match.group()
            for dict in batch_info_dicts:
                if ID in dict['Filename']:
                    newName = dict['Sample Description'][12:] + '.tiff'
                    os.rename(folder + '/' + file, folder + '/' + newName)
                    print(dict['Filename'] + ' renamed to ' + newName)

                else:
                    pass
    return

def interp_y(x0, y0, x):
    interp_func = interp1d(x0, y0, fill_value='extrapolate')

    return interp_func(x)

def bkg_subtraction(sample,bkg):
    q = sample[0]
    I= sample[1]
    s = sample[2]

    q_bkg = bkg[0]
    scale = np.average(I[-10:]/bkg[1][-10:])
    I_bkg= bkg[1]*scale
    s_bkg= bkg[2]*scale

    if q.all() != q_bkg.all():
        I_bkg = interp_y(q_bkg,I_bkg, q)
        s_bkg = interp_y(q_bkg,s_bkg, q)

    new_I = I - I_bkg
    new_s = np.abs(s - s_bkg)
    new_sample = [q, new_I, new_s]
    return new_sample

def convert_tiff_files(folder_path, x_axis_variable, binning_style, number_of_bins, bkg_path = None, show_plots = False, save_files = True, external_AgB = None, rename = False):

    if rename:
        rename_files(folder_path)


    # def lorenz_peak(xc, gam):
    def convert_AgB_array(scatter_array,center,beamstop_r=7):
        # theta, r, intensity
        masked_points = [[],[]]
        points = [[],[],[]]
        scattered = [[],[],[]]
        for y_coord, row in enumerate(scatter_array):
            for x_coord, val in enumerate(row):
                I = scatter_array[y_coord][x_coord]
                if val >= 0 and np.sqrt((x_coord-center[0]) ** 2 + (y_coord - center[1]) ** 2) > beamstop_r:
                    points[0].append(x_coord)
                    points[1].append(y_coord)
                    points[2].append(I)

                else:
                    masked_points[0].append(x_coord)
                    masked_points[1].append(y_coord)

                if val>0:
                    scattered[0].append(x_coord)
                    scattered[1].append(y_coord)
                    scattered[2].append(I)
        return points, masked_points, scattered
    def auto_mask(scatter_array, center):
        ##all of this is done in terms of pixels (cartesian pixels)

        points,masked_points,scattered = convert_AgB_array(scatter_array, center)
        def theta_mask(points,masked_points, center):
            x0 = center[0]
            y0 = center[1]

            x= np.array(points[0])
            y = np.array(points[1])
            I = np.array(points[2])
            r = np.sqrt((x-x0)**2+(y-y0)**2)
            theta = np.arctan2(y-y0,x-x0)

            theta_bins = np.linspace(-np.pi,np.pi,360)
            theta_bindex = np.digitize(theta, theta_bins)
            bin_Is = [0]*len(theta_bins)

            for i,angle in enumerate(theta):
                bin = theta_bindex[i]
                bin_Is[bin] += I[i]

            avg_I_per_angle= np.average(bin_Is)
            theta_good = []
            theta_bad = []
            r_bad = []
            good_coords = [[],[],[]]
            r_good =[]
            I_good = []
            testx=[]
            testy=[]
            beamstop = [[],[]]
            for i,angle in enumerate(theta):
                bin = theta_bindex[i]
                if bin_Is[bin]>0.015* avg_I_per_angle:
                    # bs_theta.append(angle)
                    theta_good.append(angle)
                    r_good.append(r[i])
                    I_good.append(I[i])
                    good_coords[0].append(x[i])
                    good_coords[1].append(y[i])
                    good_coords[2].append(I[i])


                else:
                    r_bad.append(np.sqrt((x[i]-x0)**2+(y[i]-y0)**2))
                    theta_bad.append(angle)
                    testx.append(x[i])
                    testy.append(y[i])
                    beamstop[0].append(x[i])
                    beamstop[1].append(y[i])

            # get beamstop width

            indx = beamstop[1] == max(beamstop[1])
            xs = [beamstop[0][i] for i,val in enumerate(indx) if val ]
            ys = [beamstop[0][i] for i, val in enumerate(indx) if val]
            bad_angles = np.array(theta_bad)
            bad_theta_range = max(bad_angles)-min(bad_angles)
            print(bad_theta_range*180/np.pi)
            print(max(r_bad))
            arclength = max(r_bad)*bad_theta_range
            ## use small angle appx




            masked_points[0] = masked_points[0]+beamstop[0]
            masked_points[1] = masked_points[1]+beamstop[1]


            return masked_points,theta_good, r_good, I_good, good_coords

        masked_points, theta_good, r_good, I_good, good_coords = theta_mask(points, masked_points, center)


        polar = [theta_good,r_good,I_good]


        return polar, masked_points, good_coords

    def radial_binning(r,I, distance, beamstop_r=7, binned_var = 'r', binning = 'linear', num_bins = 400):
        r = np.array(r)


        if binned_var == 'r':
            xaxis = "Radial pixel"
        if binned_var == 'q':
            two_theta = np.arctan(r*0.172/distance)
            r = 4*np.pi*np.sin(two_theta/2)/1.54
            beamstop_r = 4 * np.pi * np.sin(np.arctan((beamstop_r * 0.172 / distance) / 2) / 1.54)
            xaxis = 'q (1/A)'
        max_r = r.max()
        min_r = r.min()
        if binning == 'linear':
            r_bin_edges = np.linspace(beamstop_r, max_r, num = num_bins)
            r_bin_centers = (r_bin_edges[1:]+r_bin_edges[:-1])/2

        if binning == 'log':
            r_bin_edges = np.logspace(np.log10(min_r), np.log10(max_r),num = num_bins, base = 10)
            r_bin_centers = np.log10((10**r_bin_edges[1:]+10**r_bin_edges[:-1])/2)


        r_bindex = np.digitize(r, r_bin_edges)

        r_bin_I = [0]*(len(r_bin_edges))
        counts = [0]*(len(r_bin_edges))


        for i, val in enumerate(r):

            bin = r_bindex[i]-1
            r_bin_I[bin]+=I[i]
            counts[bin]+=1
        avgs=[]
        for i, total in enumerate(r_bin_I):
            if counts[i]>0:
                avgs.append(total/counts[i])
            else:
                avgs.append(0)

        ssr = [0]*(len(r_bin_edges))
        for i,val in enumerate(r):
            bin = r_bindex[i]-1
            ssr[bin]+=(I[i]-avgs[bin])**2

        ssr_array = np.array(ssr)
        count_array = np.array(counts)

        min_ppb = count_array.min()
        max_ppb = count_array.max()
        print(f'{np.average(count_array)} points per bin on average with {min_ppb} min per bin and {max_ppb} max per bin.'
              f'median: {np.median(count_array)} ')

        avgs = np.array(avgs)
        stdev_array = np.sqrt(ssr_array / count_array**2)
        ## last bin average is the amount of intensity measurements that exceed  max radius
        ## This bin is always zero by my design and binning
        avgs_for_bin_centers = (avgs[:-1])
        error_for_bin_centers = stdev_array[:-1]
        r_bin_data = [r_bin_centers, avgs_for_bin_centers, error_for_bin_centers]
        return r_bin_data
    def AgB_cal_ring_one(scatter_array):
        ## this outputs a rough center
        AgB_fig = plt.figure('Initial',figsize=(8, 8))
        ax_AgB = AgB_fig.add_subplot(111)
        ax_AgB.imshow(np.log10(scatter_array), cmap='terrain')
        ax_AgB.set_facecolor('black')
        ax_AgB.invert_yaxis()
        ring_point = plt.ginput(1)
        plt.close()

        print(ring_point)
        scatter_points = [[],[],[]]
        used_points = []
        x_ring_1=int(np.round(ring_point[0][0],0))
        y_ring_1=int(np.round(ring_point[0][1],0))

        print(x_ring_1,y_ring_1)

        used_points.append([x_ring_1, y_ring_1, scatter_array[y_ring_1][x_ring_1]])

        for i, row in enumerate(scatter_array):
            scatter_row= []
            for j, val in enumerate(row):
                if val > 0:
                    scatter_points[0].append(j)
                    scatter_points[1].append(i)
                    scatter_points[2].append(val)

        def crawl(x0,y0,scatter_grid,used_points = None):
            sys.setrecursionlimit(30000)
            I0 = scatter_grid[y0][x0]
            if used_points is None:
                used_points =[[x0,y0,I0]]


            new_cols = [x0-1,x0+1]
            new_rows = [y0-1,y0+1]

            points = []

            for row in new_rows:
                for col in new_cols:
                    I_point = scatter_grid[row][col]
                    points.append([col,row,I_point])

            viable = [point for point in points if point not in used_points and point[2]>=0.1*I0]

            if len(viable)>0:
                for point in viable:
                    used_points.append(point)
                    crawl(point[0],point[1],scatter_grid, used_points)
                else:
                    return used_points




        used_points = crawl(x_ring_1,y_ring_1,scatter_array)

        used_array = np.array(used_points)
        x_ring = used_array[:,0]
        y_ring = used_array[:,1]
        I_ring= used_array[:,2]

        sorted_indices= np.argsort(I_ring)
        x_ring1 = x_ring[sorted_indices]
        y_ring1 = y_ring[sorted_indices]
        I_ring1= I_ring[sorted_indices]

        fraction = len(x_ring)//4

        x_ring1 = x_ring1[-fraction:]
        y_ring1 = y_ring1[-fraction:]
        I_ring1 = I_ring1[-fraction:]

        points = [[x_ring1[i],y_ring1[i]] for i,val in enumerate(x_ring1)]
        fit_fig = plt.figure('AgB fitting', figsize = (8,8))
        ax_fit = fit_fig.add_subplot(111)
        ax_fit.imshow(np.log10(scatter_array),cmap='terrain')
        ax_fit.invert_yaxis()
        ax_fit.set_facecolor('black')
        for point in points:
            ax_fit.plot(point[0], point[1], ls = '', marker = 's', markersize = 1, color = 'red')

        x0,y0, r_pixel, sig =  circle_fit.hyperLSQ(points)
        r_mm = 0.172*r_pixel
        d_AgB = 58.380 * 1E-7
        # angstroms to mm
        distance = r_mm * d_AgB / (1.54 * 1E-7)





        center = [x0,y0]

        ax_fit.plot(center[0], center[1], ls = '', marker = '*', markersize = 5, color = 'red')
        plt.show()

        return center, distance, scatter_points

    def get_norm_factor(sample_name,run_dict, distance_mm, pixel_size_mm):

        for sample_dict in run_dict:

            if sample_name == sample_dict['Sample Description'][12:]:

                scalar = 1 / (pixel_size_mm / distance_mm) ** 2 / float(sample_dict['Measurement Time']) / float(sample_dict['Transmission']) / float(sample_dict['Io'])

        return scalar
    def tiff_to_polar(data_array, good_coords, scalar):
        x0 = center[0]
        y0 = center[1]


        good_data_I= []
        x = np.array(good_coords[0])
        y = np.array(good_coords[1])
        r = np.sqrt((x-x0)**2+(y-y0)**2)
        theta = np.arctan2(y-y0,x-x0)
        for i, x in enumerate(good_coords[0]):
                y = good_coords[1][i]
                good_data_I.append(data_array[y][x])
        polar_raw = [r,theta,good_data_I]
        polar_norm = [r,theta,np.array(good_data_I)*scalar]
        return polar_raw, polar_norm

    # def write_calibration_data(raw_data, norm_data, scalar, binning, binned_va)
    def write_2Dto1D(raw_data, norm_data, scalar, sample_to_detector_distance, run_list_path, folder_path, samplename,data_bkg_sub = None, bkg_name = None):

        date = datetime.datetime.now()

        dir_numpy_data = os.path.join(folder_path,"npy_files")
        dir_path_csv_data = os.path.join(folder_path,'csv_data')

        ## RAW DATA ##
        raw_array = np.array(raw_data)

        q_raw = raw_data[0]
        I_raw = raw_data[1]
        s_raw = raw_data[2]


        dir_npy_raw = os.path.join(dir_numpy_data,"raw_2Dto1D")

        if not os.path.isdir(dir_npy_raw):
            os.makedirs(dir_npy_raw)

        raw_sample_npy_path = os.path.join(dir_npy_raw,samplename+'_raw.npy')
        np.save(raw_sample_npy_path, raw_array)
        dir_path_raw = os.path.join(dir_path_csv_data, "raw_2Dto1D")

        if not os.path.isdir(dir_path_raw):
            os.makedirs(dir_path_raw)

        raw_sample_csv_path = os.path.join(dir_path_raw, samplename +'_raw.csv')

        with open(raw_sample_csv_path, 'w', newline='') as new_raw_file:
            raw_writer = csv.writer(new_raw_file, delimiter=',')
            raw_writer.writerow([f"** Data generated by WB's SAXS 2Dto1D python script {date}"])
            raw_writer.writerow([f'Sample to detector distance: {round(sample_to_detector_distance,3)} mm'])
            raw_writer.writerow(['q (1/A)', 'Intensity (arb)', 'Error (arb)'])
            for i, q in enumerate(q_raw):
                raw_writer.writerow([q, I_raw[i], s_raw[i]])

            new_raw_file.close()
        ## END RAW DATA ##





        ## NORMALIZED DATA ##
        norm_array = np.array(norm_data)
        q_norm = norm_data[0]
        I_norm = norm_data[1]
        s_norm = norm_data[2]
        dir_npy_norm = os.path.join(dir_numpy_data,"normalized_2Dto1D")
        if not os.path.isdir(dir_npy_norm):
            os.makedirs(dir_npy_norm)
        norm_sample_npy_path = os.path.join(dir_npy_norm, samplename+'_norm.npy')
        np.save(norm_sample_npy_path, norm_array)

        dir_path_norm = os.path.join(dir_path_csv_data, "normalized_2Dto1D")

        if not os.path.isdir(dir_path_norm):
            os.makedirs(dir_path_norm)

        norm_sample_csv_path = os.path.join(dir_path_norm, samplename +'_norm.csv')

        with open(norm_sample_csv_path, 'w', newline='') as new_norm_file:
            norm_writer = csv.writer(new_norm_file, delimiter=',')
            norm_writer.writerow([f"**Data generated by WB's SAXS 2Dto1D python script {date}**"])
            norm_writer.writerow([f'Data normalized by scalar value: {scalar}'])
            norm_writer.writerow([f'Scalar obtained from {run_list_path}'])
            norm_writer.writerow([f'Sample to detector distance: {round(sample_to_detector_distance,3)} mm'])
            norm_writer.writerow(['q (1/A)', 'Intensity (arb)', 'Error (arb)'])
            for i, q in enumerate(q_norm):
                norm_writer.writerow([q, I_norm[i], s_norm[i]])

            new_norm_file.close()

        if data_bkg_sub is not None:
            bkg_array = np.array(data_bkg_sub)
            q_bkg = data_bkg_sub[0]
            I_bkg = data_bkg_sub[1]
            s_bkg = data_bkg_sub[2]


            dir_npy_bkg = os.path.join(dir_numpy_data, "bkg_subtracted")
            if not os.path.isdir(dir_npy_bkg):
                os.makedirs(dir_npy_bkg)

            bkg_sample_npy_path = os.path.join(dir_npy_bkg, samplename+'_bkg_sub.npy')
            np.save(bkg_sample_npy_path,bkg_array)

            dir_path_bkg = os.path.join(dir_path_csv_data, "bkg_subtracted")
            if not os.path.isdir(dir_path_bkg):
                os.makedirs(dir_path_bkg)

            bkg_sample_csv_path = os.path.join(dir_path_bkg, samplename+'_bkg_sub.csv')

            with open(bkg_sample_csv_path, 'w', newline='') as new_bkg_csv_file:
                bkg_writer = csv.writer(new_bkg_csv_file, delimiter =',')
                bkg_writer.writerow([f"**Data generated by WB's SAXS 2Dto1D python script {date}**"])
                bkg_writer.writerow([f'Data normalized by scalar value: {scalar}'])
                bkg_writer.writerow([f'Scalar obtained from {run_list_path}'])
                bkg_writer.writerow([f'Sample to detector distance: {round(sample_to_detector_distance,3)}mm'])
                bkg_writer.writerow([f'Background subtracted using {bkg_name}'])
                bkg_writer.writerow(['q (1/A)', 'Intensity (arb)', 'Error (arb)'])

                for i, q in enumerate(q_bkg):
                    bkg_writer.writerow([q, I_bkg[i], s_bkg[i]])

                new_bkg_csv_file.close()

    x_unit_dict = {'q':'(1/A)', 'r': 'pixels', 'd': '(A)', '2θ': '(rad)'}
    units = x_unit_dict[x_axis_variable]

    files=os.listdir(folder_path)

    batch_info_dicts=[]
    data_files = []
    for file in files:
        if 'AgB' in file and file.endswith('tiff') or 'agb' in file and file.endswith('tiff') or 'Agb' in file and file.endswith('tiff'):
            print(file)
            AgB_name = file.strip('.tiff')
            AgB_file = os.path.join(folder_path, file)

        elif file.endswith('.tiff'):
            data_files.append(file)


        elif file.endswith('.csv'):
            run_list_path = os.path.join(folder_path, file)
            RunDicts=csv.DictReader(open(run_list_path))
    if external_AgB is not None:
        AgB_name = external_AgB.strip('.tiff')

        AgB_image = im.open(AgB_file)


    data_files.sort()
    for dict in RunDicts:
        batch_info_dicts.append(dict)

    if external_AgB is None:
        AgB_image = im.open(AgB_file)
        AgB_array = np.array(AgB_image)

        center, sample2detector_mm, scatter_points = AgB_cal_ring_one(AgB_array)
        polar_coords, masked_points, good_coords = auto_mask(AgB_array, center)
        scalar = get_norm_factor(AgB_name, batch_info_dicts, sample2detector_mm, 0.172)


        dir_path_centering = os.path.join(folder_path, 'AgB_centering')
        if not os.path.isdir(dir_path_centering):
            os.makedirs(dir_path_centering)


        mask_df = pd.DataFrame(np.transpose(masked_points),columns=['x','y'])
        point_df = pd.DataFrame(np.transpose(good_coords),columns=['x','y','I'])

        beam_info= os.path.join(dir_path_centering, 'beam_info' + '.txt')
        with open(beam_info, 'w') as f:
            f.write(f'Sample to detector: {sample2detector_mm}\n')
            f.write(f'center pixel coordinates: {center[0]}, {center[1]}')

        mask_df.to_csv(os.path.join(dir_path_centering, 'masked_points.csv'), index = False)
        mask_npy = mask_df.to_numpy()
        np.save(os.path.join(dir_path_centering, 'mask.npy'),mask_npy)

        point_df.to_csv(os.path.join(dir_path_centering, 'valid_pixels.csv'), index= False)
        point_npy = point_df.to_numpy()
        np.save(os.path.join(dir_path_centering, 'valid_pixels.npy'),point_npy)


        np.save(dir_path_centering, AgB_array)


        plt.savefig(os.path.join(dir_path_centering, 'Mask.png'))

        AgB_mask = plt.figure(f'mask_from_{AgB_file}')
        ax_mask = AgB_mask.add_subplot(111)
        ax_mask.invert_yaxis()
        ax_mask.set_facecolor('blue')

        ax_mask.imshow(np.log10(AgB_array), cmap = 'terrain')
        plt.plot(masked_points[0], masked_points[1], ls='', marker='s', markersize=1, color = 'black')
        plt.plot(center[0], center[1], ls='', marker='*', markersize=1, color = 'red')
    # else:
    #     AgB_import_list= os.listdir(external_AgB)
    #     for file in AgB_import_list:
    #         if file.endswith('npy')
    AgB_polar, AgB_norm_polar = tiff_to_polar(AgB_array, good_coords, scalar)

    AgB_raw = radial_binning(AgB_polar[0], AgB_polar[2], distance=sample2detector_mm, binned_var=x_axis_variable,
                             binning=binning_style, num_bins=number_of_bins)
    AgB_norm = radial_binning(AgB_norm_polar[0], AgB_norm_polar[2], distance=sample2detector_mm,
                              binned_var=x_axis_variable, binning=binning_style, num_bins=number_of_bins)

    if show_plots:
        plt.figure('raw')
        plt.errorbar(AgB_raw[0], AgB_raw[1], yerr = AgB_raw[2], ecolor ='black', label = 'AgB')
        plt.figure('norm')
        plt.errorbar(AgB_norm[0], AgB_norm[1], yerr=AgB_norm[2], ecolor='black', label = 'AgB')

    scatter_path = os.path.join(folder_path, 'log10_scatter_plots')
    if not os.path.isdir(scatter_path):
        os.makedirs(scatter_path)
    for file in data_files:
        name = file.removesuffix('.tiff')
        data_path = os.path.join(folder_path, file)
        data_image = im.open(data_path)
        data_array = np.array(data_image)
        log_img = im.fromarray(np.log10(data_array))
        log_img.save(os.path.join(scatter_path,file))
        scalar = get_norm_factor(name, batch_info_dicts, sample2detector_mm, 0.172)
        AgB_polar,AgB_norm_polar = tiff_to_polar(data_array, good_coords, scalar)
        data_norm = radial_binning(AgB_norm_polar[0], AgB_norm_polar[2], distance=sample2detector_mm,
                                   binned_var=x_axis_variable, binning=binning_style, num_bins=number_of_bins)
        data_raw = radial_binning(AgB_polar[0], AgB_polar[2], distance=sample2detector_mm, binned_var=x_axis_variable,
                                  binning=binning_style, num_bins=number_of_bins)
        if bkg_path is None and save_files:
            write_2Dto1D(data_raw, data_norm, scalar, sample2detector_mm, run_list_path,folder_path,name)
        elif save_files:
            bkg_name = os.path.abspath(bkg_path)
            bkg = np.load(bkg_path)
            data_bkg_sub = bkg_subtraction(data_norm, bkg)
            write_2Dto1D(data_raw, data_norm, scalar, sample2detector_mm, run_list_path, folder_path, name, data_bkg_sub, bkg_name = bkg_name )



        if show_plots:
            fig_raw = plt.figure('raw',x_scale = 'log')
            plt.errorbar(data_raw[0], data_raw[1], yerr = data_raw[2], ecolor ='black', label = name)

            fig_norm = plt.figure('norm')
            plt.errorbar(data_norm[0], data_norm[1], yerr = data_norm[2], ecolor ='black', label = name)

            if bkg_path is not None:
                fig_bkg = plt.figure('background_subtracted')
                plt.errorbar(data_bkg_sub[0], data_bkg_sub[1], yerr=data_bkg_sub[2], ecolor = 'black', label = name, ls = '', marker = '.', alpha = 0.5)


    if show_plots:

        plt.figure('raw')
        plt.xlabel(x_axis_variable+' '+units)
        plt.ylabel("I (arb)")
        plt.xscale('log')
        plt.yscale('log')
        fig_raw.legend(loc="outside upper right")


        plt.figure('norm')
        plt.xlabel(x_axis_variable+' '+units)
        plt.ylabel("I (arb)")
        plt.xscale('log')
        plt.yscale('log')
        fig_norm.legend(loc="outside upper right")

        if bkg_path is not None:
            plt.figure('norm')
            plt.errorbar(bkg[0],bkg[1],yerr=bkg[2],ecolor = 'black',label = 'bkg')

            plt.figure('background_subtracted')
            plt.errorbar(bkg[0], bkg[1], yerr=bkg[2], ecolor='black', label='bkg')
            plt.xlabel(x_axis_variable+' '+units)
            plt.ylabel("I (arb)")
            plt.xscale('log')
            plt.yscale('log')
            fig_bkg.legend(loc="outside upper right")


        plt.show()


if __name__ == '__main__':
#### TEST CASE ####
    folder = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20250626_aldy3_size_gel_stocks'
    bkg_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/Backgrounds/SAXS_Hexane_bkg_norm.npy'
    # AgB_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20241120_WB1176_synth/SAXS_AgB.tiff'
    convert_tiff_files(folder, 'q','log',400, bkg_path=bkg_path, save_files = True, show_plots=True, external_AgB=None, rename = False)

