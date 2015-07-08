#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np, matplotlib as mp, scipy as sp
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from functions import comp, fits

# READ FILE

file_name = sys.argv[1]

data = np.genfromtxt(file_name, dtype="float", unpack = True)

z = data[0,:]
rich_log = np.log10(data[1,:])
mass = data[2,:]

# BIN

min_z = min(z)
max_z = max(z)
z_bin_size = 0.5
n_z_bins = int((max_z - min_z) / z_bin_size)

min_rich = min(rich_log)
max_rich = max(rich_log)
rich_bin_size = 0.5
n_rich_bins = int((max_rich - min_rich) / rich_bin_size)

# DEFAULTS

min_mass = 13.0
max_mass = 15.0
mass_bin_size = 0.5
n_mass_bins = int((max_mass - min_mass) / mass_bin_size)

# EXTRAC XY

def data_limit(data, bin_size, i):
    limit_low = min(data) + i * bin_size
    limit_up = limit_low + bin_size
    limit = ((data >= limit_low) & (data < limit_up))
    return limit

def get_err(x, y):
    xerr = 1.0 / (np.sqrt(x) * np.log(10))
    yerr = np.zeros(len(y)) + 0.01
    return xerr, yerr

def extract_xy(x_data, y_data, z_data, y_bin, z_bin, i, j):
    limit = ((data_limit(y_data, y_bin, j)) & (data_limit(z_data, z_bin, i)))
    x = x_data[limit]
    y = y_data[limit]
    xerr, yerr = get_err(x, y)
    return x, y

# PLOT

fig = plt.figure()
gs = gridspec.GridSpec(n_z_bins, n_rich_bins)
fig.subplots_adjust(hspace = 0.0001, wspace = 0.0001)

k = -1

ax = []
for i in range(n_rich_bins * n_z_bins):
    ax.append(0)

for i in range(n_z_bins - 1, -1, -1):
    for j in range(n_rich_bins):
        k += 1
        if j == 0:
            ax[k] = fig.add_subplot(gs[k])
        else:
            ax[k] = fig.add_subplot(gs[k])
            ax[k].yaxis.set_visible(False)
        if i != n_z_bins - 1:
            ax[k].xaxis.set_visible(False)
            
        # EXTRACT DATA
        rich_limit_low = min_rich + j * rich_bin_size
        rich_limit_up = rich_limit_low + rich_bin_size
        rich_limit = ((rich_log >= rich_limit_low) & (rich_log < rich_limit_up))
        z_limit_low = min_z + i * z_bin_size
        z_limit_up = z_limit_low + z_bin_size
        z_limit = ((z >= z_limit_low) & (z < z_limit_up))
        limit = ((rich_limit) & (z_limit))

        #if len(rich_log[limit]) > 0:
         #   y = (rich_log[limit] - min(rich_log[limit])) / (max(rich_log[limit]) - min(rich_log[limit]))
        #else:
        y1 = rich_log[limit]
        
        x1 = mass[limit]

        #x, y = extract_xy(mass, rich_log, z, rich_bin_size, z_bin_size, i, j)

        #print y
        #print y1
        
        # PLOT DATA
        ax[k].scatter(x1, y1, s=2, c='0', lw=0)
        ax[k].set_xlim(min_mass, max_mass)
        ax[k].set_ylim(0.0, 1.0)
        ax[k].xaxis.set_major_locator(MaxNLocator(3))
        ax[k].yaxis.set_major_locator(MaxNLocator(3))
fig_name = file_name + '_test_plot.1.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name
