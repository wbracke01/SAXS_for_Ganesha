import numpy as np
import pandas as pd
from PIL import Image as im
import matplotlib.pyplot as plt
import matplotlib



path = '/Users/willbrackett/Library/CloudStorage/OneDrive-TheUniversityofTexasatAustin/Research/SAXS_data/20251006_PAPEG_20nm_B2_series/WB2020_17nm_PAPEG_240mgmL.tiff'

image = im.open(path)

data = np.log10(np.array(image))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', ['blue','green', 'yellow', 'red'])
fig, ax = plt.subplots()
ax.set_facecolor("darkblue")         # Set axes background color
ax.imshow(data, cmap = cmap)

plt.show()