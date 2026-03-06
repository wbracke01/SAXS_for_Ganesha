from scipy.interpolate import interp1d
from scipy.optimize import curve_fit as cf
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import os
import csv
import pandas as pd
from matplotlib.widgets import SpanSelector





def csv_to_data_array(sample_path):
    if sample_path.endswith('.csv'):
        sample_header = []
        with open(sample_path,'r') as f:
            csv_reader = csv.reader(f)
            end_header=False
            q = []
            I = []
            s = []

            for i, row in enumerate(csv_reader):
                if row[0] != 'q (1/A)':
                    if end_header:
                        q.append(float(row[0]))
                        I.append(float(row[1]))
                        s.append(float(row[2]))
                    else:
                        sample_header.append(row[0])
                else:
                    end_header = True

            q = np.array(q)
            I = np.array(I)
            s = np.array(s)
            sample_array = [q,I,s]


        return sample_array, sample_header


def interp_y(x0, y0, x):
    interp_func = interp1d(x0, y0, fill_value='extrapolate')

    return interp_func(x)
def find_first_peak(sample, name):
    er = sample[2]
    x = sample[0]
    y = sample[1]
    max_y = max(y)

    x = x[er<0.01*max_y]
    y = y[er<0.01*max_y]
    ## doing this in kind of a dirty linear way even though these points are log space bc i only care about the log values
    def box_averaging(x,y, half_box):
        xavg = []
        yavg = []
        for i, val in enumerate(x[half_box:-half_box]):
            xavg.append(np.mean(x[i-half_box:i+half_box]))
            yavg.append(np.mean(y[i-half_box:i+half_box]))
        xavg = np.array(xavg)
        yavg = np.array(yavg)
        return (xavg,yavg)

    xavg, yavg = box_averaging(x,y,10)
    xlog = np.log10(xavg)
    ylog = np.log10(yavg)


    x1= np.log10((10**xlog[1:]+10**xlog[:-1])/2)
    y1 = np.log10((10**ylog[1:]+10**ylog[-1:])/2)
    d1 = np.diff(ylog)/np.diff(xlog)*x1

    d1_clean = d1[~np.isnan(d1)]
    x1_clean = x1[~np.isnan(d1)]

    min_d1_clean = np.min(d1_clean)

    crit_x = 10**x1_clean[d1_clean == min_d1_clean]
    fig,ax = plt.subplots(2,1, sharex = True)
    fig.suptitle(name)
    ax[0].set_title('Minimum q for fitting')
    ax[0].set_xlabel('q (1/A)')
    ax[0].set_ylabel('I (arb)')
    ax[0].plot(xavg,yavg)
    ax[0].vlines(crit_x, ymax = max(y), ymin = min(y), color = 'red')
    ax[0].set_yscale('log')

    ax[1].set_title('Differential representation')
    ax[1].set_xlabel('q (1/A)')
    ax[1].set_ylabel('dI/dq')
    ax[1].plot(10**x1,10**d1)
    ax[1].plot(crit_x, 10**min_d1_clean, ls = '', marker = '^', markersize = 5, color = 'red')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    fig.tight_layout()



    return crit_x


def struct(sample, form, q_min, q_max):
    ## ABSOLUTELY MAKE SURE YOU USE BACKGROUND SUBTRACTED DATA




    q = sample[0]
    I = sample[1]
    s = sample[2]

    q_f = form[0]
    I_f = form[1]



    I_f = interp_y(q_f, I_f, q)


    ## values denoted as hi should be q values that exist past the first peak


    q_hi = np.array([val for val in q if val>q_min and val<q_max])
    I_hi = np.array([I[i] for i,val in enumerate(q) if val>q_min and val<q_max])


    I_f_hi = np.array([I_f[i] for i,val in enumerate(q) if val>q_min and val<q_max])


    scaled_hi_q = I_f_hi / I_hi




    noninf = scaled_hi_q[scaled_hi_q<np.inf]

    struct_scaling = np.average(noninf)

    print(struct_scaling)

    I_scaled = struct_scaling*I
    s_scaled = struct_scaling*s

    I_struct = I_scaled/I_f



    sample_struct = [q,I_struct]
    sample_scaled = [q,I_scaled]



    return sample_struct, sample_scaled

def single_struct(sample_path,form_path, q_range = None):
    if form_path.endswith('csv'):
        form_array, form_header = csv_to_data_array(form_path)
    elif form_path.endswith('npy'):
        form_array = np.load(form_path)
    else:
        raise ValueError('Form factor should be a file from 2Dto1D code')

    if sample_path.endswith('csv'):
        sample_array, sample_header = csv_to_data_array(sample_path)
    elif sample_path.endswith('npy'):
        sample_array = np.load(sample_path)
    else:
        raise ValueError("Sample data should be files processed by Will's 2Dto1D code")

    name = sample_path.split()[-1]



def batch_struct_factors(folder, form_path, q_range = None):
    if form_path.endswith('csv'):
        form_array, form_header = csv_to_data_array(form_path)
    elif form_path.endswith('npy'):
        form_array = np.load(form_path)
    else:
        raise ValueError('Form factor should be a file from 2Dto1D code')

    q = form_array[0]
    I = form_array[1]
    # I_error = form_array[:,2]
    fig1, ax1 = plt.subplots(subplot_kw={'xscale': 'log', 'yscale': 'log'})
    ax1.plot(q, I, marker='o', markerfacecolor='none', color='dodgerblue', ls='', zorder=0, label='data')

    if q_range is None:
        def onselect(vmin, vmax):
            print(vmin, vmax)
            plt.close()

        span = SpanSelector(ax1, onselect, 'horizontal',
                            props=dict(facecolor='blue', alpha=0.5))
        plt.show()
        x1, x2 = span.extents


    fig_struct = plt.figure('struct')
    ax_struct = fig_struct.add_subplot(111)
    fig_fit = plt.figure('overlay')
    ax_fit = fig_fit.add_subplot(111)
    ax_fit.plot(form_array[0],form_array[1], ls = '', marker = 'o', markersize = 4, color = 'black', fillstyle = 'none')
    try:
        files = os.listdir(folder)

        files.sort()
    except:
        files = [folder]


    for file in files:
        if not file.endswith('.DS_Store'):
            file_path = os.path.join(folder, file)
            print(file)
            name = os.path.basename(file)
            print(name)
            name = name.strip('bkg_sub')
            print(name)

            if file.endswith('csv'):
                sample_array, sample_header = csv_to_data_array(file_path)
            elif file.endswith('npy'):
                sample_array = np.load(file_path)
            else:
                raise ValueError("Sample data should be files processed by Will's 2Dto1D code")
            sample_struct, sample_scaled = struct(sample_array,form_array, x1,x2)
            ax_struct.plot(sample_struct[0], sample_struct[1], ls = '',
                               marker = 'o', markersize = '4', label = name, alpha = 0.5)
            ax_fit.plot(sample_scaled[0],sample_scaled[1],label = name)

    ax_struct.legend(loc = 'upper right')

    fig_fit.legend(loc = 'outside upper right')
    ax_struct.set_xlabel('q (1/A)')
    ax_struct.set_ylabel('S(q)')
    ax_struct.set_xlim(q[0], q[-1])
    ax_struct.set_ylim(-10,10)
    fig_raw, ax_raw = plt.subplots()
    ax_raw.plot(sample_array[0], sample_array[1], ls='', marker='o', label=name)
    ax_raw.plot(q,I, marker='o', ls = '', label = 'form')
    fig_raw.legend(loc='outside upper right')
    fig_raw.suptitle('raw')

    df_sample_scaled = pd.DataFrame({'q': sample_scaled[0], 'I_scaled': sample_scaled[1]})
    df_sample_scaled.to_csv('sample_scaled.csv', index=False)

    # Save form_array
    df_form_array = pd.DataFrame({'q_f': form_array[0], 'I_f': form_array[1]})
    df_form_array.to_csv('form_array.csv', index=False)

    plt.show()
    return



sample_pathnpy = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/B2_charge_size_series/10nm/8p_10,6nm_c2_100.npy'
form_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/B2_charge_size_series/10nm/8p_10,6nm_1mgmL.npy'
single = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20250618_11_28_gel_stocks/npy/bkg_sub/11nm_gel_stock.npy'
batch_struct_factors(sample_pathnpy, form_path)
