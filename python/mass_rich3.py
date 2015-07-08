#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np, matplotlib as mp, scipy as sp
from scipy import stats
from scipy.odr import *
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

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
        rich_limit_low = min_rich + j * rich_bin_size
        rich_limit_up = rich_limit_low + rich_bin_size
        rich_limit = ((rich_log >= rich_limit_low) & (rich_log < rich_limit_up))
        z_limit_low = min_z + i * z_bin_size
        z_limit_up = z_limit_low + z_bin_size
        z_limit = ((z >= z_limit_low) & (z < z_limit_up))
        limit = ((rich_limit) & (z_limit))
        x = rich_log[limit]
        y = mass[limit]
        ax[k].scatter(x, y, s=2, c='0', lw=0)
        ax[k].set_xlim(min_rich, max_rich)
        ax[k].set_ylim(min_mass, max_mass)
        ax[k].xaxis.set_major_locator(MaxNLocator(3))
        ax[k].yaxis.set_major_locator(MaxNLocator(3))

#ax.set_ylabel('Log$_{10}$ M$_{Halo}$')
#ax.set_xlabel('Log$_{10}$ N$_{gal}$')
#ax.set_title(file_name + " mass-richness relation")
fig_name = file_name + '_mass_v_rich_plot.1.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name
