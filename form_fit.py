import csv

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import pandas as pd
import matplotlib

from scipy.optimize import curve_fit as cf
from scipy.interpolate import interp1d
import os
import periodictable as pt 
from periodictable import xsf


import pretty_errors


def ito_contrast(tinfrac,bkg = None):


    bd = {'dmf':['C3H7NO',.944], 'hexane': ['C6H8',.661]}
    

    indium = 1-tinfrac
    tin = tinfrac
    ox = 3/2*indium+2*tin

    total = np.sum([indium,tin,ox])

    in_frac = indium/total
    sn_frac = tin/total
    o_frac = ox/total

    unit = f'In{in_frac}Sn{sn_frac}O{o_frac}'
    
    
   

    compound = pt.formula('In2O3')

    compound.density = 7.18

    nc_sld, mu = xsf.xray_sld(unit, density = 7.14, energy = 8.04)  #E10 / cm^2

    if bkg is not None:
        bkg_sld, mu = xsf.xray_sld(bd[bkg][0], density = bd[bkg][1], energy = 8.04) #E10/cm^2

        print(nc_sld-bkg_sld)
        return((nc_sld - bkg_sld)**2)

def oz_fit(q,s0,z):
    return(s0/(1-z*q**2))

def quadratic(q,s0,n):
    return s0+n*q**2


def init_r_fit_number(q,r,scale):
    # lengths/V
    V = 4 / 3 * np.pi * r ** 3
    F = 3*(np.sin(q * r) - q * r * np.cos(q * r)) / (q * r) ** 3
    I = scale*np.abs(F)**2*V**2
    return I

def form_fit_number(q, r_mean, sig, scale):
    q = np.atleast_1d(q)
    # radius grid
    rmax = r_mean + 3.5 * sig
    rmin = r_mean - 3.5 * sig
    r_bins = np.linspace(rmin, rmax, 1000)
    r_bins = 0.5 * (r_bins[1:] + r_bins[:-1])
    delta_r = r_bins[1] - r_bins[0]

    # number-weighted Gaussian
    D_n = (1 / (sig * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((r_bins - r_mean) / sig) ** 2)

    # broadcasted q*r matrix (N_r, N_q)
    r_col = r_bins[:, None]
    q_row = q[None, :]
    qr = r_col * q_row

    # form factor, handle qr -> 0 limit (F -> 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        F = 3 * (np.sin(qr) - qr * np.cos(qr)) / (qr ** 3)
    F = np.where(np.isfinite(F), F, 1.0)

    # volumes (N_r,1) and assemble intensity contributions
    V = (4.0 / 3.0) * np.pi * (r_col ** 3)
    I_r = (np.abs(F) ** 2) * (V ** 2) * (D_n[:, None] * delta_r)

    # sum over radii and scale
    summation = scale * np.sum(I_r, axis=0)
    return summation


    # scale is lengths/Φ
    rmax = r_mean + 3.5 * sig
    rmin = r_mean - 3.5 * sig
    r_bins = np.linspace(rmin, rmax, 1000)
    r_bins = 0.5 * (r_bins[1:] + r_bins[:-1])
    delta_r = r_bins[1] - r_bins[0]

    # Number-weighted Gaussian distribution
    D_n = (1 / (sig * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((r_bins - r_mean) / sig) ** 2)

    V_r = (4/3)*np.pi*r_bins**3

    # Volume-weighted distribution: D_V = V × D_N / ∫ V D_N
    D_v = V_r * D_n
    D_v /= np.sum(D_v * delta_r)  # normalize
    terms=[]
    for i, r in enumerate(r_bins):
        qr = q*r
        V = 4 / 3 * np.pi * r ** 3

        F = 3*((np.sin(qr)-qr*np.cos(qr))/qr**3)

        I = V*F ** 2 * D_v[i] * delta_r
        terms.append(I)
    term_array = np.array(terms)
    summation = scale * np.sum(term_array, axis=0)
    return summation

def get_sample(path):
    if path.endswith('.csv'):
        data = pd.read_csv(path, header = 0)
        data = data.dropna()
        q = np.array(data["q (1/A)"])
        I = np.array(data["Intensity (arb)"])
        I_error = np.array(data["Error (arb)"])
    if path.endswith('.npy'):
        data = np.load(path)
        q = data[0,:]
        I = data[1,:]
        I_error = data[2,:]

    good_indx = np.where((np.isfinite(I)) & (np.isfinite(I_error)))

    q = q[good_indx]
    I = I[good_indx]
    I_error = I_error[good_indx]

    return q, I, I_error

def update_form_log(log, fit_list):

    def sort_log_by_diameter(log):
        with open(log, 'r') as f:
            lines = f.readlines()
            def extract_diameter(line):
                try:
                    for part in line.split('|'):
                        if 'D_mean' in part:
                            return float(part.split(':')[1].split()[0])
                except Exception:
                    return float('inf')  # Put malformed lines at the end
            f.close()
            sorted_lines = sorted(lines, key=extract_diameter)

            with open(log, 'w') as f:
                f.writelines(sorted_lines)
                f.close()

    

    log_path = log
    for fit in fit_list:
        if fit['fit_params'][1]/fit['fit_params'][0]>0.1:
            name = fit['name'] + '(x)'
        else:
            name = fit['name']
        print(name)
        new_line = (
            f"{name:30s} | "
            f"D_mean: {2*fit['fit_params'][0]/10:8.2f} nm | "
            f"D_stdev: {2*fit['fit_params'][1]/10:8.2f} nm | "
            f"PDI: {fit['fit_params'][1]/fit['fit_params'][0]:8.2f} |"
            f"Scale: {fit['fit_params'][2]:.3e} | "
            f"VolFrac: {fit['volfrac']*100:6.3f}% | "
            f"Conc: {fit['mass_conc']:8.2f} mg/mL"
        )
        if 'mass_conc' in list(fit.keys()):
            new_line += f" | volfrac: {fit['volfrac']:8.2f}"
            new_line += f" | mass_conc: {fit['mass_conc']:8.2f} mg/mL"

        # Read existing lines (if any)
        existing_lines = []
        if os.path.exists(log_path):
            with open(log_path, 'r') as f:
                existing_lines = f.readlines()

    # Check for existing entry and update or append
        updated = False
        for i, line in enumerate(existing_lines):
            if line.startswith(name):
                existing_lines[i] = new_line + "\n"
                updated = True
                break

        if not updated:
            existing_lines.append(new_line + "\n")

    # Write all lines back to file
        with open(log_path, 'w') as f:
            f.writelines(existing_lines)

            f.close()

    sort_log_by_diameter(log)

    print('information logged')

def conc_est_with_form(f, tin_frac, bkg, cap_thickness = 0.088):

    ito_density = 7140 # mg/mL

    rho_squared = ito_contrast(tin_frac, bkg)*1E20 # cm^-4
    r_a = f['fit_params'][0] #A
    r_cm = r_a*(1E-10/1E-2)
    v_p = 4/3*np.pi*r_cm**3 # cm^3

    cap_thickness = cap_thickness # cm

    q_gun = np.logspace(-5,-3,100)
    I_gun = form_fit_number(q_gun,*f['fit_params']) #adjust for thickness in cm

    def guinier(q,I0):
        return I0 * np.exp(-q ** 2 * r_a ** 2 / 3)
    
    p, cov = cf(guinier,q_gun,I_gun)
    I0 = p[0]/cap_thickness
  
    num_density = I0/(v_p*rho_squared)
    volfrac = num_density * v_p

    mass_conc = ito_density * num_density

    return volfrac, mass_conc
    

def run_single_fit(path,distribution = 'number'):

    path_end = path.split('/')[-1]
    samplename = path_end.split('.')[0]
    print(samplename)

    if distribution =='number':
        init = init_r_fit_number
        final = form_fit_number
    else:
        raise ValueError('Distribution must be "number" or "volume"')
    
    q,I,I_error = get_sample(path)

    r_guess = 70
    sig_guess = 7.3
    scale_guess = 1E-5
    rmax = r_guess*1.5
    rmin = r_guess*0.5
    sigmax =r_guess
    sigmin = 0.01*r_guess

    pi, ci = cf(init, q, I, p0 = [r_guess,scale_guess], bounds =((0,0),(np.inf,np.inf)))
    r_init = pi[0]
    scale_init = pi[1]

    fig1, ax1 = plt.subplots(subplot_kw={'xscale':'log', 'yscale':'log'})
    fig1.suptitle(samplename + ' form fitting')
    ax1.errorbar(q,I, I_error, marker = 'o', markerfacecolor = 'none', color = 'dodgerblue', ls = '', zorder = 0, label ='data')

    def onselect(vmin, vmax):
        plt.close()

    span = SpanSelector(ax1, onselect, 'horizontal',
                                 props=dict(facecolor='blue', alpha=0.5))
    plt.show()
    x1,x2 = span.extents

    fit_q = q[np.where((q>x1)&(q<x2))]
    fit_I = I[np.where((q>x1)&(q<x2))]
    fit_I_error = I_error[np.where((q>x1)&(q<x2))]

    p, c = cf(final, fit_q, fit_I, maxfev=10000, p0 = [r_init, 0.1* r_init, scale_init], bounds = ((0,0,0), (np.inf,np.inf,np.inf)), sigma = fit_I_error)

    I_fitted = final(q, *p)

    fit_dict = {'name': samplename, 'fitting': distribution, 'q_I_form':[q,I,I_fitted], 'fit_params': p, 'fit_span': (x1, x2)}

    return fit_dict


def run_form_fitting(data, distribution, save=None, plot = True, estimate_conc=False, tin_frac = 0.05, solvent = 'hexane' ):
    def plot_form_fits(fits):

        if len(fits)>1:
            fig, ax = plt.subplots(subplot_kw={'xscale': 'log', 'yscale': 'log'})
            ax.set_xlabel('q (1/A)')
            ax.set_ylabel('Intensity (arb)')
        
        if not isinstance(fits,list):
            fits = [fits]

        plot_colors = mcolors.LinearSegmentedColormap.from_list('mycmap', ['purple', 'dodgerblue'], N=len(fits))


        for i,f in enumerate(fits):
            fig_single, ax_single = plt.subplots(subplot_kw={'xscale': 'log', 'yscale': 'log'})
            q = f['q_I_form'][0]
            I_sample = f['q_I_form'][1]
            I_fitted = f['q_I_form'][2]

            ax_single.errorbar(q, I_sample, marker='o', markerfacecolor='none', color='dodgerblue', ls='', zorder=0, label='data')
            ax_single.plot(q, I_fitted, color='red', label='fitted')
            ax_single.vlines(f['fit_span'],min(I_sample),max(I_sample), ls='--', color='red', label = 'fitting range')

            ax_single.set_xlabel('q (1/A)')
            ax_single.set_ylabel('Intensity (arb)')
            ax_single.text(min(q), min(I_fitted), f'{f["fitting"]} distribution\nD: {2*f["fit_params"][0] / 10:.4} nm\nStDev: {2*f["fit_params"][1] / 10:.4} nm\nScale: {f["fit_params"][2]:.2e}')

            ax_single.legend()
            ax_single.set_title(f['name']+' form')

            if len(fits)>1:
                ax.errorbar(q, I_sample, marker='o', markerfacecolor='none', color=plot_colors(i), ls='', zorder=0, label=f['name'])
                ax.plot(q, I_fitted, color='red', ls='--')

            if save is not None:
                fig_path = os.path.join(save, 'figures', 'form_factor')

                if not os.path.exists(fig_path):
                    os.makedirs(fig_path)

                savename = os.path.join(fig_path, f['name'] + '_form_fit.svg')
                fig_single.savefig(savename)
        if len(fits)>1:
            
            ax.legend()

        plt.show()

    def save_form_factor(fits): 

        for f in fits:
            print(f['name'])
            csv_path = os.path.join(save,'form_factor', 'csv')
            npy_path = os.path.join(save,'form_factor', 'npy')
            for loc in [csv_path, npy_path]:
                if not os.path.exists(loc):
                    os.makedirs(loc)

            csv_file = os.path.join(csv_path, f['name'] + '_form_fit.csv')
            npy_file = os.path.join(npy_path, f['name'] + '_form_fit.npy')

            # Log the form factor data
            with open(csv_file, 'w') as cf:
                cf.write('q (1/A) ,I_sample (arb), I_form (arb) \n')
                for q_val, I_sample, I_form in zip(f['q_I_form'][0], f['q_I_form'][1], f['q_I_form'][2]):
                    cf.write(f'{q_val},{I_sample},{I_form}\n')

            np.save(npy_file, f['q_I_form'])

    data_dicts = []

    if os.path.isdir(data):
        print('running batch')

        for filename in os.listdir(data):

            
            if filename.endswith('.npy'):
                path = os.path.join(data, filename)
                samplename = path.split('/')[-1]
                savename = samplename.split('.')[0]

                fit_dict = run_single_fit(path, distribution)
                data_dicts.append(fit_dict)

    
    else:
        if data.endswith('.npy'):
            path = data
            samplename = path.split('/')[-1]
            savename = samplename.split('.')[0]

            fit_dict = run_single_fit(path, distribution)
            data_dicts.append(fit_dict)
    
    if estimate_conc:
        for f in data_dicts:
            f['volfrac'], f['mass_conc'] = conc_est_with_form(f, tin_frac, solvent)

            print(f['name'], f['volfrac'], f['mass_conc'])

    if save is not None:
        save_form_factor(data_dicts)

    if plot:
        plot_form_fits(data_dicts)


    return data_dicts

def structure_factor(data, save=None,same_particles = False, estimate_conc=False, tin_frac = 0.08, solvent = 'dmf', show_plots = False):

    s_dicts = []
    

    if os.path.isdir(data):
        data_files = os.listdir(data)

        colors = mcolors.LinearSegmentedColormap.from_list('mycmap', ['orange', 'blue'], N=len(data_files))

        for filename in data_files:
            ds = {}

            if filename.endswith('.npy'):
                ds['file_path'] = os.path.join(data, filename)
                path = os.path.join(data, filename)
                samplepath = path.split('/')[-1]
                name = samplepath.split('.')[0]

                ds['q'], ds['I'], ds['I_error'] = get_sample(path)
                ds['avg_I'] = np.mean(ds['I'])

                ds['name'] = name

                s_dicts.append(ds)
    else:
        ds = {}
        ds['file_path'] = data
        path = data
        samplepath = path.split('/')[-1]
        name = samplepath.split('.')[0]

        ds['q'], ds['I'], ds['I_error'] = get_sample(path)
        ds['avg_I'] = np.mean(ds['I'])

        ds['name'] = name

        s_dicts.append(ds)
    
    def individual_particles_structure_factor(s_dicts):

        for f in s_dicts:
            f['form_fit'] = run_single_fit(f['file_path'], distribution = 'number')
        
        fig_struct, ax_struct = plt.subplots(subplot_kw={'xscale': 'log'})
        ax_struct.set_xlabel('q (1/A)')
        ax_struct.set_ylabel('S(q)')
        for f in s_dicts:
            q = f['q']
            I_sample = f['I']
            form = f['form_fit']['q_I_form'][2]

            f['s'] = I_sample/form
            f['d'] = f['form_fit']['fit_params'][0]*2
            if f['form_fit']:
                f['volfrac'], f['mass_conc'] = conc_est_with_form(f['form_fit'], tin_frac, solvent)
                label = f['name'] + f' {f["mass_conc"]:.2f} mg/mL'
            else:
                label = f['name']

            ax_struct.plot(q, f['s'], marker = 'o', ls = '', label = label)
        ax_struct.legend()
        ax_struct.set_ylim(-5,5)


        if show_plots:
            plt.show()

                
   
    def same_particle_structure_factor(s_dicts):
        
        cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', ['orange','green', 'blue'], N = len(s_dicts))
        if same_particles: 
            print('using common form factor from lowest concentration')

            s_dicts = sorted(s_dicts, key=lambda x: x['avg_I'])


            lowest_conc_fit = run_single_fit(s_dicts[0]['file_path'], distribution = 'volume')
            lowest_c_q = lowest_conc_fit['q_I_form'][0]
            lowest_c_form = lowest_conc_fit['q_I_form'][2]

            if estimate_conc:
                lowest_v, lowest_c = conc_est_with_form(lowest_conc_fit, tin_frac, solvent)


            fig_scale, ax_scale = plt.subplots(subplot_kw={'xscale': 'log', 'yscale': 'log'})
            ax_form = ax_scale.twinx()
            ax_form.set_yscale('log')
            ax_form.set_ylabel('I form (arb)')
            ax_form.tick_params(axis='y', colors = 'red')

            ax_form.plot(lowest_c_q, lowest_conc_fit['q_I_form'][1], color = 'red', ls = '', marker = 'o')
            ax_form.plot(lowest_c_q, lowest_c_form, color = 'black', ls = '--', label = 'form factor')
        
            ax_scale.set_xlabel('q (1/A)')
            ax_scale.set_ylabel('I/I form')

            for i, f in enumerate(s_dicts):

                name = f['name']
                q = f['q']
                I_sample = f['I']

                interp_form = interp1d(lowest_c_q, lowest_c_form, fill_value = 'extrapolate')

                ax_scale.plot(q, I_sample/interp_form(q), color = cmap(i), marker = 'o', ls = '', label = name)
            
            def onclick(qmin,qmax):
                plt.close(fig_scale)

            q_span = SpanSelector(ax_scale, onclick, 'horizontal')
            plt.show()

            for i, f in enumerate(s_dicts):
                q = f['q']
                I_sample = f['I']
                q1,q2 = q_span.extents
                f['span'] = (q1,q2)

                q_in_range = q[(q >= q1) & (q <= q2)]
                I_sample_scalar = I_sample[(q >= q1) & (q <= q2)]

                I_common_interp = interp1d(lowest_c_q, lowest_c_form, fill_value = 'extrapolate')

                scalar = np.mean(I_sample_scalar/I_common_interp(q_in_range))
                f['scalar'] = scalar
                f['s'] = I_sample/(scalar*I_common_interp(q))
                f['d'] = lowest_conc_fit['fit_params'][0]*2

                if estimate_conc:
                    f['volfrac'] = lowest_v*scalar
                    f['mass_conc'] = lowest_c*scalar
                    particle_volume_nm = 4/3*(f['d']/20)**3 * np.pi

                    f['num_dens'] = f['volfrac'] / particle_volume_nm

            
            fig_struct, ax_struct = plt.subplots(subplot_kw={'xscale': 'log'})
            ax_struct.set_xlabel('q (1/A)')
            ax_struct.set_ylabel('S(q)')

            fig_structqd, ax_structqd = plt.subplots(subplot_kw={'xscale': 'log'})
            ax_structqd.set_xlabel('qd (unitless)')
            ax_structqd.set_ylabel('S(qd)')

            fig_real, ax_real = plt.subplots(subplot_kw={'xscale': 'log'})
            
            for i, f in enumerate(s_dicts):
                label = f['name']

                q = f['q']
                I_sample = f['I']
                scalar = f['scalar']
                d = f['d']
                qd = q*d
                s = f['s']
                real = 2*np.pi/(q)/10

                if estimate_conc: 
                    label =  f'{f['mass_conc']:.0f} mgmL'
                    title = f'{f['d']/10:.1f} nm ' + label
                else:
                    label = f['name']
                    title = f['name']


                ax_struct.plot(q, s, color = cmap(i), marker = 'o', ls = '', label = label)
                ax_structqd.plot(qd, s, color = cmap(i), marker = 'o', ls = '', label = label)
                ax_real.plot(real, s, color = cmap(i), marker = 'o', ls = '', label = label)

                fig_s1, ax_s1 = plt.subplots(subplot_kw={'xscale': 'log'})
                ax_s1.set_xlabel('q (1/A)')
                ax_s1.set_ylabel('S(q)')
                ax_s1.set_title(title)

                roi = (qd>0)&(qd<9.1)
                q_roi = q[roi]
                s_roi = s[roi]

                ax_s1.plot(q, s, color = cmap(i), marker = 'o', ls = '', label = label)

                ax_s1.set_xlim(min(q_roi), max(q_roi))
                ax_s1.set_ylim(min(s_roi)*0.9, max(s_roi)*1.1)

                fig_qd1, ax_qd1 = plt.subplots(subplot_kw={'xscale': 'linear'})
                ax_qd1.set_xlabel('qd (unitless)')
                ax_qd1.set_ylabel('S(qd)')
                ax_qd1.set_title(title)

                roi_qd = (qd>0)&(qd<9.1)    
                qd_roi = qd[roi_qd]
                s_roi = s[roi_qd]
        
                ax_qd1.plot(qd, s, color = cmap(i), marker = 'o', ls = '', label = label)

                ax_qd1.set_xlim(min(qd_roi), max(qd_roi))
                ax_qd1.set_ylim(min(s_roi)*0.9, max(s_roi)*1.1)

                fig_real1, ax_real1 = plt.subplots(subplot_kw={'xscale': 'log'})
                ax_real1.set_xlabel(r'$\frac{2 \pi}{q} (nm)$')
                ax_real1.set_ylabel(r'$S(\frac{2 \pi}{q})$')
                ax_real1.set_title(title)

                ax_real1.plot(real, s, color = cmap(i), marker = 'o', ls = '', label = label)
                ax_real1.set_xlim(min(real[roi_qd]), max(real[roi_qd]))
                ax_real1.set_ylim(min(s_roi)*0.9, max(s_roi)*1.1)



                if save is not None:

                    new_paths = ['q','qd','real']

                    save_paths = []
                    for p in new_paths:
                        fig_path = os.path.join(save, 'figures', 'structure_factors',p)
                        save_paths.append(fig_path)

                        if not os.path.exists(fig_path):
                            os.makedirs(fig_path)

                    savename = os.path.join(save_paths[0], f['name'] + '_s(q).svg')
                    fig_s1.savefig(savename)
                    savename = os.path.join(save_paths[1], f['name'] + '_s(r).svg')
                    fig_qd1.savefig(savename)
                    savename = os.path.join(save_paths[2], f['name'] + '_real.svg')
                    fig_real1.savefig(savename)
            
            ax_struct.legend()
            ax_structqd.legend()
            ax_real.legend()

            if save is not None:
                fig_path = os.path.join(save, 'figures', 'structure_factors')

                if not os.path.exists(fig_path):
                    os.makedirs(fig_path)

                savename = os.path.join(fig_path, 'all_s(q).svg')
                fig_struct.savefig(savename)
                savename = os.path.join(fig_path, 'all_s(qd).svg')
                fig_structqd.savefig(savename)

        if show_plots:
            plt.show()

        plt.close('all')
        return s_dicts, lowest_conc_fit

    if same_particles:
        return same_particle_structure_factor(s_dicts)

    else:
        return individual_particles_structure_factor(s_dicts)

def s0_fit_with_q_limit(q,s, start_q):

    index = q<start_q
    q = q[index]
    s = s[index]

    pi, ci = cf(quadratic, q, s)
    fiti = quadratic(q,*pi)

    below_1 = np.where(fiti < 0.9)
    q_new = q[below_1]
    s_new = s[below_1]
    
    if len(q_new) < 2:
        return pi, ci, q
    else: 
        p_new, c_new = cf(quadratic, q_new, s_new, maxfev = 1000)

    return p_new, c_new, q_new

def b2_plotting(b2_dicts):


    fig, ax = plt.subplots(subplot_kw={'xscale': 'log'})
    ax.set_ylabel('S(q)')
    ax.set_xlabel('q (1/A)')
    colors = mcolors.LinearSegmentedColormap.from_list('mycmap', ['orange', 'green', 'blue'], N=len(b2_dicts))
    

    for i, db2 in enumerate(b2_dicts):
        scalar = db2['scalar']

        ax.plot(db2['q'][:-1], db2['s'][:-1], marker='o', alpha = 0.5, color = colors(i), ls='', label = f'{db2['mass_conc']:.2f}')
    
    for i, db2 in enumerate(b2_dicts):
        ax.plot(db2['q_s0'], quadratic(db2['q_s0'], db2['s0'], db2['z']), color = colors(i))

    q = db2['q']

    roi = np.where((q*db2['d']>0)&(q*db2['d']<9))

    ax.set_xlim(min(q[roi]), max(q[roi]))
    ax.set_ylim(0.5,1.75)

    ax.legend()

    return fig

def b2_same_particles(b2_dicts, form_dilute, save = None, save_structure = False):

    

    colors = mcolors.LinearSegmentedColormap.from_list('mycmap', ['orange', 'green', 'blue'], N=len(b2_dicts))

    r_g_a = form_dilute['fit_params'][0]  # A
    r_g_nm = r_g_a / 10 # nm
    b2_hs = 16/3*np.pi*r_g_nm**3 # nm^3


    def get_x_from_click(ax=None):
        if ax is None:
            ax = plt.gca()
        x_selected = {}
        def onclick(event):
            if event.inaxes == ax:
                x_selected['x'] = event.xdata
                plt.close()
        cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        
        return x_selected.get('x')

    
    for db2 in b2_dicts:
        
    
        fig_max_q, ax_max_q = plt.subplots(subplot_kw={'xscale': 'log'})
        ax_max_q.plot(db2['q'][:-1], db2['s'][:-1], marker='o', markerfacecolor='none', ls='', color = 'black')
        ax_max_q.set_ylim(0,1.2)
        x_val = get_x_from_click(ax_max_q)
        db2['max_q'] = x_val
        print("Selected x:", x_val)

        q_max = x_val
        index = db2['q'] < q_max
        q_new = db2['q'][index]
        s0_new = db2['s'][index]
        
        s0_fit, s0_cov= cf(quadratic, q_new, s0_new)
        
        db2['q_s0'] = q_new
        print(q_new)
        

        db2['s0'] = s0_fit[0]
        db2['s0_err'] = np.sqrt(np.diag(s0_cov))[0]

        db2['z'] = s0_fit[1]
        db2['z_err'] = np.sqrt(np.diag(s0_cov))[1]

    struct_fit_fig = b2_plotting(b2_dicts)

    fig_b2, ax_b2 = plt.subplots(1,2, layout = 'tight')
    fig_b2.suptitle('$S(q) \\approx S(0) + nq^2$')
    ax_b2[0].set_title(f'{2/10*form_dilute['fit_params'][0]:.2f} nm')
    ax_b2[0].set_xlabel('$\\rho_n (1/nm)$')
    ax_b2[0].set_ylabel('$\\frac{1}{S(0)}-1$')

    ax_b2[1].set_title(f'{2/10*form_dilute['fit_params'][0]:.2f} nm')
    ax_b2[1].set_xlabel('$\\rho_n (1/nm)$')
    ax_b2[1].set_ylabel('n')

    fig_qclick, ax_qclick = plt.subplots(figsize=(5,5), layout='tight')
    ax_qclick.set_ylabel('q clicked position (1/A)')
    ax_qclick.set_xlabel('Sample concentration (1/nm)')

    fit_c = []
    fit_s_inv = []
    fit_num_dens = []

    for i, db2 in enumerate(b2_dicts):
        if db2['s0_err']<db2['s0']:
            fit_num_dens.append(db2['num_dens'])
            fit_s_inv.append(1/db2['s0']-1)
            s0 = db2['s0']
            s0_err = db2['s0_err']
            ax_b2[0].errorbar(db2['num_dens'],1/s0-1,1/s0*s0_err/s0, marker='o', color =colors(i), ls='')
            ax_b2[1].errorbar(db2['num_dens'],np.abs(db2['z']), db2['z_err'], marker='o', color =colors(i), ls='')

            ax_qclick.plot(db2['num_dens'], db2['max_q'], marker='o', color =colors(i), ls='', label = f'{db2['mass_conc']:.2f} mg/mL')

    def s0_conc_line(x,a):
        return 2*a*x
    def s0_conc_line_int(x,b,c):
        return c+2*b*x
    def s0_conc_power(x,d,alpha):
        return d*x**alpha


    def r_squared(y_true, y_pred):
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - ss_res / ss_tot

    s_inv_array = np.array(fit_s_inv)

    # change fitting range basedon behavior
    p_line, cov_line= cf(s0_conc_line, fit_num_dens[1:], s_inv_array[1:])

    p_line_int, cov_line_int = cf(s0_conc_line_int, fit_num_dens, s_inv_array)

    p_power, cov_power = cf(s0_conc_power, fit_num_dens, s_inv_array)



    smooth_conc = np.linspace(0, max(fit_num_dens), 1000)
    dens_array = np.array(fit_num_dens)

    fit_line = s0_conc_line(dens_array, *p_line)
    fit_line_int = s0_conc_line_int(dens_array, *p_line_int)
    fit_line_power = s0_conc_power(dens_array, *p_power)


    r2_line = r_squared(fit_s_inv, fit_line)
    r2_line_int = r_squared(fit_s_inv, fit_line_int)
    r2_line_power = r_squared(fit_s_inv, fit_line_power)

    a = p_line[0]
    a_err = np.sqrt(np.diag(cov_line))[0]
    b = p_line_int[0]
    b_err = np.sqrt(np.diag(cov_line_int))[0]
    c = p_line_int[1]
    c_err = np.sqrt(np.diag(cov_line_int))[1]
    d = p_power[0]
    d_err = np.sqrt(np.diag(cov_power))[0]
    alpha = p_power[1]
    alpha_err = np.sqrt(np.diag(cov_power))[1]

    ax_b2[0].plot(smooth_conc, s0_conc_line(smooth_conc, *p_line), color='red', ls='--', label='linear')
    # ax_b2[0].plot(smooth_conc, s0_conc_line_int(smooth_conc, *p_line_int), color='green', ls='--', label='linear (intercept)')
    # ax_b2[0].plot(smooth_conc, s0_conc_power(smooth_conc, *p_power), color='blue', ls='--', label='power')

    6/3*np.pi*r_g_nm**3

    fit_text = (
        f'Linear: 1/S(0)-1 = 2a ρₙ\n'
        f'  a = {a:.3g} ± {a_err:.2g}\n'
        f'  B₂/B₂ₕₛ = {a/b2_hs:.3f} ± {a_err/b2_hs:.3f}\n'
        f'  $d_eff$ = {(a*3/2/np.pi)**(1/3):.2f} nm\n'
        f'  $R^2$ = {r2_line:.2f}\n'
    )
    ax_b2[0].text(
        0.05, 0.95,
        fit_text,
        transform=ax_b2[0].transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='left',
        bbox=dict(boxstyle='round', facecolor='none', edgecolor='none', alpha=0.7)
    )
    ax_b2[0].legend()
    # ax_qclick.legend()
    if save:
        b2_path = os.path.join(save, 'figures', 'b2_fits')
        if not os.path.exists(b2_path):
            os.makedirs(b2_path, exist_ok=True)
        
        struct_fit_path = os.path.join(b2_path, 'structure_fits.svg')
        struct_fit_fig.savefig(struct_fit_path)
        b2_fit_path = os.path.join(b2_path, 'b2_fits.svg')
        fig_b2.savefig(b2_fit_path)
        qclick_path = os.path.join(b2_path, 'q_clicks.svg')
        fig_qclick.savefig(qclick_path) 


def save_fit(fit_list, save_path):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    fit_file = os.path.join(save_path, 'form_fits.npy')

    np.save(fit_file, fit_list)
    print('fits saved')

if __name__ == '__main__':

    matplotlib.style.use('/Users/brackw/Library/CloudStorage/Dropbox-UniversityofMichigan/Will Brackett/Python_projects/Readable_style.txt')
    # folder = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20251003_PAPEG_10nm_B2_series/npy/bkg_sub'
    # save_title = 'PAPEG_10'

    data_path= '/Users/willbrackett/Downloads/SAXS Data/npy/bkg_sub'
    save_title = None
    log_path = 'number_average_sheet.txt'

    fit_list = run_form_fitting(data_path, 'number', estimate_conc = True, tin_frac=0.05, solvent='dmf', save = save_title)
    update_form_log(log_path, fit_list)

    s_dicts = structure_factor(data_path, save=save_title, same_particles=False, estimate_conc=True, tin_frac=0.05, solvent='hexane', show_plots=True)
    # b2_same_particles(s_dicts, form_dilute, save = save_title)


