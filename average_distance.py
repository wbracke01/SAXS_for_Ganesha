import numpy as np 
import matplotlib.pyplot as plt
import matplotlib



def number_density(d_nm,c_mg_per_mL):
    ITO_density = 7140 #mg/mL

    volume_particle_nm3 = (4/3)*np.pi*(d_nm/2)**3 #nm^3
    volume_particle_mL = volume_particle_nm3 * 1e-21 #mL
    mass_particle_mg = ITO_density * volume_particle_mL #mg
    n_particles_per_mL = c_mg_per_mL / mass_particle_mg
    n_particles_per_nm3 = n_particles_per_mL * 1e-21 #particles/nm^3
    
    return n_particles_per_nm3

def average_spherical_distance(d_nm,c_mg_per_mL):
    rho_n = number_density(d_nm,c_mg_per_mL)

    r_avg = (4/3 *np.pi*rho_n)**(-1/3)
    return r_avg

if __name__ == "__main__":

    matplotlib.rcParams.update({'font.size': 16})
    d = [10,15,20]

    fig, ax = plt.subplots(1,2, layout = 'tight')
    fig_q, ax_q = plt.subplots()

    for d in d:
        c = np.linspace(0,200,100) #mg/mL
        r_avg = average_spherical_distance(d,c) #nm
        rho_n = number_density(d,c) 

        ax[0].plot(c,rho_n,label=f'd={d} nm')
        ax[1].plot(c,r_avg,label=f'd={d} nm')

        ax_q.plot(c, 2*np.pi/(r_avg*10), label=f'd={d} nm')

    

    ax[0].set_xlabel('Concentration (mg/mL)')
    ax[0].set_ylabel('Number Density (particles/nm^3)')
    ax[0].legend()
    ax[1].set_xlabel('Concentration (mg/mL)')
    ax[1].set_ylabel('Average Spherical Distance (nm)')
    ax[1].legend()

    ax_q.set_xlabel('Concentration (mg/mL)')
    ax_q.set_ylabel('q corresponding to max avg distance (1/A)')
    ax_q.legend()


    plt.show()