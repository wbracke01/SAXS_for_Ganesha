import numpy as np
import matplotlib.pyplot as plt
import csv
import pretty_errors
from matplotlib.widgets import SpanSelector
import periodictable as pt

def read_csv_data(file_path):
    header = {}
    data = {}
    q = []
    I = []

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0].startswith('#'):
                header[row[1]] = row[2]  # Store header info

            else:    
                try:
                    q.append(float(row[0]))
                    I.append(float(row[1]))
                except (ValueError, IndexError):
                    continue  # Skip rows with invalid data
    data['q'] = np.array(q)
    data['I'] = np.array(I)
    return data, header

def SiO2_contrast(tinfrac,bkg = None):

    
    
   

    compound = pt.formula('In2O3')

    compound.density = 0.14

    nc_sld, mu = xsf.xray_sld(unit, density = 7.14, energy = 8.04)  #E10 / cm^2

    if bkg is not None:
        bkg_sld, mu = xsf.xray_sld(bd[bkg][0], density = bd[bkg][1], energy = 8.04) #E10/cm^2

        return((nc_sld - bkg_sld)**2)


def fit_log_linear(q, I):
    log_q = np.log10(q)
    log_I = np.log10(I)

    coeffs = np.polyfit(log_q, log_I, 1)
    slope, intercept = coeffs

    return slope, intercept

def select_multiple_spans(data):
    q, I = data['q'], data['I']
    spans = []

    fig, ax = plt.subplots()
    ax.plot(q, I)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('q (1/Å)')
    ax.set_ylabel('I')

    def onselect(qmin, qmax):
        spans.append((qmin, qmax))
        ax.axvspan(qmin, qmax, color='blue', alpha=0.3)
        fig.canvas.draw()

    span = SpanSelector(ax, onselect, 'horizontal', props=dict(facecolor='blue', alpha=0.3))

    print("Select spans. Close the window when done.")

    spans.sort(key=lambda x: x[0])  # Sort spans by starting q-value
    plt.show()
    return spans

def fit_linear_regions(data):
    q = data['q']
    I = data['I']
    
    spans = select_multiple_spans(data)
    print(spans)
    indxs = [np.where((q >= span[0]) & (q <= span[1])) for span in spans]
    q_ranges = [q[indx] for indx in indxs]
    I_ranges = [I[indx] for indx in indxs]

    slopes = [slope for slope, intercept in [fit_log_linear(q_range, I_range) for q_range, I_range in zip(q_ranges, I_ranges)]]

    fig, ax = plt.subplots()
    ax.set_xlabel('q (1/Å)')
    ax.set_ylabel('I')

    ax.plot(q,I, color = 'black')
    for i, (q_range, I_range, slope) in enumerate(zip(q_ranges, I_ranges, slopes)):
        ax.plot(q_range, I_range, linewidth = 3, label = f'slope: {slope:.2f}')

    ax.legend()
    plt.show()


def plot_raw(data):
    fig_Iq, ax_Iq = plt.subplots(figsize=(6, 4))

    q = data['q']
    I = data['I']
    ax_Iq.plot(q, I, marker='o', linestyle='None', markersize=3, color='blue')
    ax_Iq.set_xscale('log')
    ax_Iq.set_yscale('log')
    ax_Iq.set_xlabel('q (1/Å)')
    ax_Iq.set_ylabel('I')

    fig_porod, ax_porod = plt.subplots(figsize=(6, 4))
    ax_porod.plot(q, I*q**4, marker='o', linestyle='None')

    # dIdq = np.gradient(np.log(I), np.log(q))
    # fig_dIdq, ax_dIdq = plt.subplots(figsize=(6, 4))
    # ax_dIdq.plot(np.log(q), dIdq, ls = 'none', marker = '.', color='red')
    # ax_dIdq.set_ylim(-10,10)

if __name__ == '__main__':
    file_path = '/Users/willbrackett/Downloads/Andres Final Results/SAXS plot/Pure_silica_samples_B_0_00001.CSV'  # Replace with your CSV file path

    data, header = read_csv_data(file_path)
    # fit_linear_regions(data)

    fit_linear_regions(data)
    plt.show()