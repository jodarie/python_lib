#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np, matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from functions import comp, fits

# READ FILE

file_name = sys.argv[1]

data = np.genfromtxt(file_name, dtype = "float", unpack = True)

z = data[0,:]
mass = data[1,:]
rich_log = np.log10(data[2,:])

# BIN

min_z = min(z)
max_z = max(z)
z_bin_size = 0.2
n_z_bins = comp.num_bins(min_z, max_z, z_bin_size)

min_rich = min(rich_log)
max_rich = max(rich_log)
rich_bin_size = 0.2
n_rich_bins = comp.num_bins(min_rich, max_rich, rich_bin_size)

print 'min rich:', min_rich
print 'max rich:', max_rich
print 'n bins:', n_rich_bins

# DEFAULTS

min_mass = 13.0
max_mass = 15.0

# EXTRAC XY

def data_limit(data, bin_size, i):
    limit_low = min(data) + i * bin_size
    limit_up = limit_low + bin_size
    limit = ((data >= limit_low) & (data < limit_up))
    return limit

def extract_xy(x_data, y_data, z_data, x_bin, z_bin, i, j):
    limit = ((data_limit(x_data, x_bin, j)) & (data_limit(z_data, z_bin, i)))
    return x_data[limit], y_data[limit]

def scale(data):
    data = np.array(data, dtype = 'float')
    return (data - min(data)) / (max(data) - min(data))

def hist_xy(data, n_bins, min_val, max_val):
    y, x = np.histogram(data, bins = n_bins, range = (min_val, max_val))
    x = x + ((x[-1] - x[0]) / (len(x) - 1)) / 2
    x = x[np.where(x < x[-1])]
    return x, scale(y)
    
# PLOT

mp.rcParams.update({'font.size': 6})

fig = plt.figure()
fig.set_size_inches(8.27, 11.69)
gs = gridspec.GridSpec(n_z_bins, n_rich_bins)
fig.subplots_adjust(hspace = 0.0001, wspace = 0.0001)

k = -1
n_hist_bins = 10

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
        if i > 0:
            ax[k].xaxis.set_visible(False)
        x, y = extract_xy(rich_log, mass, z, rich_bin_size, z_bin_size, i, j)
        x, y = hist_xy(y, n_hist_bins, min_mass, max_mass)    
        if j == n_rich_bins - 1:
            text = 'z = ' + str(round(min_z + z_bin_size * (i + 0.5), 2))
            ax[k].plot(x, y, label = text)
            ax[k].legend()
        else:
            ax[k].plot(x, y)
        ax[k].set_xlim(min_mass, max_mass)
        ax[k].xaxis.set_major_locator(MaxNLocator(3))
        ax[k].yaxis.set_major_locator(MaxNLocator(3))
        
plt.figtext(0.47, 0.06, "Mass", fontdict = {'fontsize':18})

fig_name = file_name + '_hist_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

# PRINT DATA

data_array = []

for i in range(n_z_bins):
    for j in range(n_rich_bins):
        x, y = extract_xy(rich_log, mass, z, rich_bin_size, z_bin_size, i, j)
        x, y = hist_xy(y, n_hist_bins, min_mass, max_mass)
        if i == 0 and j == 0:
            data_array.append(x)
        data_array.append(y)

data_array = np.transpose(data_array)
a, b = data_array.shape

out_file_name = file_name + '_bin_data.txt'                                        
out_file = open(out_file_name,'w')

for i in range(a):
    for j in range(b):
        print >> out_file, data_array[i, j],
    print >> out_file, ""
print "Data saved to:", out_file_name
