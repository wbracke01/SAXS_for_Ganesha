
from matplotlib import pyplot as plt
import numpy as np
import circle_fit
from PIL import Image as im
import sys
from scipy.ndimage import center_of_mass


def agb_cal_ring_one(scatter_array):
    ## this outputs a rough center
    AgB_fig = plt.figure('Initial', figsize=(8, 8))
    ax_AgB = AgB_fig.add_subplot(111)
    ax_AgB.imshow(np.log10(scatter_array), cmap='terrain')
    ax_AgB.set_facecolor('black')
    ax_AgB.invert_yaxis()
    ring_point = plt.ginput(1)
    plt.close()

    print(ring_point)
    scatter_points = [[], [], []]
    used_points = []
    x_ring_1 = int(np.round(ring_point[0][0], 0))
    y_ring_1 = int(np.round(ring_point[0][1], 0))

    print(x_ring_1, y_ring_1)

    used_points.append([x_ring_1, y_ring_1, scatter_array[y_ring_1][x_ring_1]])

    for i, row in enumerate(scatter_array):
        scatter_row = []
        for j, val in enumerate(row):
            if val > 0:
                scatter_points[0].append(j)
                scatter_points[1].append(i)
                scatter_points[2].append(val)

    def crawl(x0, y0, scatter_grid, used_points=None):
        sys.setrecursionlimit(30000)
        I0 = scatter_grid[y0][x0]
        if used_points is None:
            used_points = [[x0, y0, I0]]

        new_cols = [x0 - 1, x0 + 1]
        new_rows = [y0 - 1, y0 + 1]

        points = []

        for row in new_rows:
            for col in new_cols:
                I_point = scatter_grid[row][col]
                points.append([col, row, I_point])

        viable = [point for point in points if point not in used_points and point[2] >= 0.1 * I0]

        if len(viable) > 0:
            for point in viable:
                used_points.append(point)
                crawl(point[0], point[1], scatter_grid, used_points)
            else:
                return used_points

    used_points = crawl(x_ring_1, y_ring_1, scatter_array)

    used_array = np.array(used_points)
    x_ring = used_array[:, 0]
    y_ring = used_array[:, 1]
    I_ring = used_array[:, 2]

    sorted_indices = np.argsort(I_ring)
    x_ring1 = x_ring[sorted_indices]
    y_ring1 = y_ring[sorted_indices]
    I_ring1 = I_ring[sorted_indices]

    fraction = len(x_ring) // 4

    x_ring1 = x_ring1[-fraction:]
    y_ring1 = y_ring1[-fraction:]
    I_ring1 = I_ring1[-fraction:]

    points = [[x_ring1[i], y_ring1[i]] for i, val in enumerate(x_ring1)]
    fit_fig = plt.figure('AgB fitting', figsize=(8, 8))
    ax_fit = fit_fig.add_subplot(111)
    ax_fit.imshow(np.log10(scatter_array), cmap='terrain')
    ax_fit.invert_yaxis()
    ax_fit.set_facecolor('black')
    for point in points:
        ax_fit.plot(point[0], point[1], ls='', marker='s', markersize=1, color='red')

    x0, y0, r_pixel, sig = circle_fit.hyperLSQ(points)

    r_mm = 0.172 * r_pixel
    d_AgB = 58.380 * 1E-7
    # angstroms to mm
    distance = r_mm * d_AgB / (1.54 * 1E-7)

    center = [x0, y0]

    ax_fit.plot(center[0], center[1], ls='', marker='*', markersize=5, color='red')
    plt.show()

    return center, distance, scatter_points


agb_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/Backgrounds/SAXS_Toluene_90mins/SAXS_AgB.tiff'
bkg_path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/Backgrounds/SAXS_Toluene_90mins/SAXS_AgB.tiff'
agb = im.open(agb_path)
agb_a = np.array(agb)

# center, distance, scatter_points = agb_cal_ring_one(agb_a)
# print(center)
center = [267.4778276232713,358.5887571060092]

row_num, col_num = agb_a.shape

Y=np.linspace(0,row_num-1,row_num)
X=np.linspace(0,col_num-1,col_num)

y,x = np.meshgrid(Y,X,indexing='ij')
x_new = x-center[0]
y_new = y-center[1]
r = np.sqrt(x_new**2+y**2)
theta = np.arctan2(y_new,x_new)
grid = np.array([agb_a, x, y,r,theta])

x_sum = np.sum(agb_a,axis=0)
y_sum = np.sum(agb_a,axis=1)

bin_list = np.linspace(0,np.max(r),400)

digitized_bin_list = np.digitize(r,bin_list)
sums = [0]*len(bin_list)

for i,bin in enumerate(bin_list):
    good = np.where(digitized_bin_list==i)
    sums[i] = np.sum(agb_a[good])

plt.plot(bin_list, sums)
plt.show()
