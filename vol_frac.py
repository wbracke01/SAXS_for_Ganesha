import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector


def form_fit_volume_fraction(path, rho, d):
    data = np.load(path)
    q = data[0,:]
    I = data[1,:]
    e = data[2,:]

    print(q)
    def onclick(minq,maxq):
        plt.close()
        return



    fig, ax = plt.subplots()

    span = SpanSelector(ax, onclick, 'horizontal',)
    ax.plot(q,I)
    plt.show()




    print(data.shape)

if __name__ == "__main__":
    delt_rho_sqr = 1420.9 #Delta Rho Squared [10^20 сm^-4]
    diameter = 10.4 #nm
    thickness = 0.098 #cm
    path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20250617_size_gel_form_factors/npy/bkg_sub/27nm_gel_stock_5mgmL.npy'
    form_fit_volume_fraction(path, delt_rho_sqr, 10.4)