import pandas as pd
from PIL import Image as im
import numpy as np

from matplotlib import pyplot as plt
import os
import csv
import pandas as pd

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


    F = 3 * (np.sin(qr) - qr * np.cos(qr)) / (qr ** 3)

    # volumes (N_r,1) and assemble intensity contributions
    V = (4.0 / 3.0) * np.pi * (r_col ** 3)
    I_r = (np.abs(F) ** 2) * (V ** 2) * (D_n[:, None] * delta_r)

    # sum over radii and scale
    summation = scale * np.sum(I_r, axis=0)
    return summation




if __name__ == '__main__':
    r_mean = 50 #A
    r_list = np.linspace(r_mean - 0.1*r_mean,r_mean+0.1*r_mean,100)
    q = np.linspace(1e-3,3e-1,1000)
    

    sig_list = 0.1
    d_rho_sqrd = 1420.9 #10^20 cm^-4

    scale = 1e-8

    I = form_fit_number(q, 50,sig_list*50,scale)
    plt.plot(q,I, label = f's = {sig_list}', marker = '.')
    plt.xscale('log')
    plt.yscale('log')

    df = pd.DataFrame({'q': q, 'I': I})
    df.to_csv('form_factor_for_anna.csv', index=False)
    plt.show()