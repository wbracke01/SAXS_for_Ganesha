import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import pandas as pd
import matplotlib


def get_data(file_path):
    pandas_df = pd.read_csv(file_path,delimiter = '\t', comment='#', header=0)

    q = pandas_df.iloc[:,0].values
    I = pandas_df.iloc[:,1].values
    sigma = pandas_df.iloc[:,2].values
    
    return q, I, sigma

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

if __name__ == '__main__':
    folder = '/Users/brackw/Library/CloudStorage/Dropbox-UniversityofMichigan/Will Brackett/Research/ligand_stripped_aps_USAXS'
    files = os.listdir(folder)
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["blue", "red"])
    colors = cmap(np.linspace(0, 1, len(files)))
    for i, file in enumerate(sorted(files)):
        print(file)
        if file.endswith('.dat'):
            file_path = os.path.join(folder, file)
            q, I, sigma = get_data(file_path)
            ax.plot(q, I,label = file.split('.')[0], color=colors[i])
    
    ax.set_xlabel('q')
    ax.set_ylabel('I')
    ax.legend()
    plt.show()