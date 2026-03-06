from form_test import ito_contrast, run_single_fit, form_fit_volume, conc_est_with_form

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def conc_est_with_form(f, tin_frac, bkg):
    ito_density = 7140 # mg/mL

    rho_squared = ito_contrast(tin_frac, bkg)*1E20 # cm^-4
    r_a = f['fit_params'][0] #A
    r_cm = r_a*(1E-10/1E-2)
    v_p = 4/3*np.pi*r_cm**3 # cm^3

    cap_thickness = 0.098 # cm

    q_gun = np.logspace(-5,-3,100)
    I_gun = form_fit_volume(q_gun,*f['fit_params']) #adjust for thickness in cm

    def guinier(q,I0):
        return I0 * np.exp(-q ** 2 * r_a ** 2 / 3)
    
    p, cov = cf(guinier,q_gun,I_gun)
    I0 = p[0]/cap_thickness
  
    volfrac = I0/(v_p*rho_squared)
    mass_conc = ito_density * volfrac

    return volfrac, mass_conc
    

def b2_same_particles(data, form ):
    b2_dicts = []

    if os.path.isdir(data):
        print('running batch')

        for filename in os.listdir(data):
            db2 = {}

            if filename.endswith('.npy'):
                db2['file_path'] = os.path.join(data, filename)
                path = os.path.join(data, filename)
                samplepath = path.split('/')[-1]
                name = samplepath.split('.')[0]

                db2['q'], db2['I'], db2['I_error'] = get_sample(path)
                db2['avg_I'] = np.mean(db2['I'])

                db2['name'] = name

               

                b2_dicts.append(db2)

    ## sort by appx concentration

    b2_dicts.sort(key=lambda x: x['avg_I'])

    b2_dicts[0]['scalar'] = 1.0

    form_dilute = run_single_fit(b2_dicts[0]['file_path'])

    form_vol_frac, form_mass_conc =conc_est_with_form(form_dilute, tin_frac = 0.08, bkg = 'dmf')

    plt.plot(form_dilute['q_I_form'][0], form_dilute['q_I_form'][1], marker='o', markerfacecolor='none', color='dodgerblue', ls='', label='data')
    plt.plot(form_dilute['q_I_form'][0], form_dilute['q_I_form'][2], color='black', ls='--', label='fit')

    print(form_vol_frac, form_mass_conc)
    b2_dicts[0]['volfrac'] = form_vol_frac
    b2_dicts[0]['mass_conc'] = form_mass_conc
    fig, ax = plt.subplots(subplot_kw={'xscale': 'log', 'yscale': 'log'})

    for db2 in b2_dicts[1:]:
        form_I = form_fit_volume(db2['q'], *form_dilute['fit_params'])
        scalar = np.average(db2['I'] / form_I, weights = db2['I_error']**-2)
        db2['scalar'] = scalar
        db2['volfrac'] = form_vol_frac * db2['scalar']
        db2['mass_conc'] = form_mass_conc * db2['scalar']
        ax.plot(db2['q'], db2['I']/ (scalar*form_I))


if __name__ == '__main__':
    folder = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20250819_B2_lig_stripped_10/npy/bkg_sub'

    