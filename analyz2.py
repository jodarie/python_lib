#! /Users/Rowen/Documents/Library/anaconda/bin/python

###########
# ANALYZ2 #
###########

import sys, argparse, os.path
import numpy as np
from scipy import stats as ss
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

# READ ARGUMENTS

parser = argparse.ArgumentParser('ANALYZ2 OPTIONS:')
parser.add_argument('-i', '--input_file', dest = 'input_file',
                    help = 'Input file name.')
parser.add_argument('-o', '--output_file', dest = 'output_file',
                    help = 'Output file name.')
parser.add_argument('-s', '--zspec_col', dest = 'zspec_col',
                  type = int, help = 'Spectroscopic redshift column number.')
parser.add_argument('-p', '--zphot_col', dest = 'zphot_col',
                  type = int, help = 'Photometric redshift column number.')
parser.add_argument('-b', '--n_bins', dest = 'n_bins', default = [200, 200],
                  nargs = '+', type = int, help = 'Number of bins' + \
                  '[Default = 200 200]')
parser.add_argument('-x', '--x_range', dest = 'x_range', type = float,
                  nargs = '+', help = 'X-axes range. (z_spec)')
parser.add_argument('-y', '--y_range', dest = 'y_range', type = float,
                  nargs = '+', help = 'Y-axes range. (z_phot)')
parser.add_argument('-y2', '--y2_range', dest = 'y2_range', type = float,
                  nargs = '+', help = 'Y-axes range. (z_phot_err)')
parser.add_argument('--mean', dest = 'mean', action = 'store_true',
                    help = 'Plot mean values.')
parser.add_argument('--fit', dest = 'best_fit', action = 'store_true',
                    help = 'Plot best fit line.')
parser.add_argument('--diff', dest = 'diff', action = 'store_true',
                    help = 'Plot z_spec - z_phot.')
parser.add_argument('-e', '--err_col', dest = 'zphoterr_col', type = int,
                    help = 'Photometric redshift error column number.')
parser.add_argument('--interp', dest = 'interp', default = 'none',
                    help = 'Bin interpolation type. [Default = none]')
opts = parser.parse_args()

if not opts.input_file:
    print 'Error: Input file not specified.'
    exit()
 
if not os.path.exists(opts.input_file):
    print 'Error: Input file not found.'
    exit()

if opts.zspec_col < 1 or opts.zphot_col < 1:
    print 'Error: Column number must be > 0.'
    exit()

if opts.n_bins[0] < 1:
    print 'Error: The number of bins must be > 0.'
    exit()

# READ FILE

data = np.genfromtxt(opts.input_file, unpack = True,
                     dtype = 'S')
z_spec = np.array(data[opts.zspec_col - 1, :], dtype = 'float')
z_phot = np.array(data[opts.zphot_col - 1, :], dtype = 'float')
if opts.zphoterr_col:
    z_phot_err = np.array(data[opts.zphoterr_col - 1, :], dtype = 'float')

# RMS

def get_scatter(x, y, z):
    df = (x - y) / (1 + x)
    df2 = np.sum(df[np.abs(df) > z] ** 2)
    return np.sqrt(df2 / len(x))
    
full_scatter = get_scatter(z_spec, z_phot, 0.0)
ol_scatter = get_scatter(z_spec, z_phot, 0.15)
print 'Full Scatter =', full_scatter
print 'Scatter minus outliers = ', ol_scatter

# START PLOT Z_SPEC VS Z_PHOT

if opts.diff:
    hist1, xedges, yedges = np.histogram2d(z_spec,
                                           (z_spec - z_phot),
                                           bins = opts.n_bins)
else:
    hist1, xedges, yedges = np.histogram2d(z_spec, z_phot,
                                           bins = opts.n_bins)

plt.figure(1)
extent1 = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
hist1[hist1 == 0.0] = None
plt.imshow(hist1.transpose()[::-1], extent = extent1,
           interpolation = opts.interp, cmap = cm.jet,
           aspect = 'auto', norm = LogNorm())
plt.colorbar()
plt.autoscale(False)
if not opts.x_range:
    plt.xlim(min(z_spec), max(z_spec))
    x = np.arange(min(z_spec), max(z_spec) + 0.1, 0.05)
else:
    plt.xlim(opts.x_range)
    x = np.arange(opts.x_range[0], opts.x_range[1] + 0.1, 0.05)
if not opts.y_range:
    if not opts.diff:
        plt.ylim(min(z_phot), max(z_phot))
    else:
        plt.ylim(-0.1, 0.1)
else:
    plt.ylim(opts.y_range)

# LINE FIT

if opts.best_fit:

    def line_func(x, m, c):
        return m * x + c

    m, c, r_value, p_value, std_err = ss.linregress(z_spec, z_phot)

    print 'Slope: %0.3f' % m
    print 'Intercept: %0.3f' % c
    print 'R Value: %0.3f' % r_value
    print 'P Value: %0.3f' % p_value
    print 'Standard Error: %0.3f' % std_err

    label = r'Best Fit: $\mu$ = %0.3f' % m + '; $b$ = %0.3f' % c
    
    plt.plot(x, line_func(x, m, c), color = 'red', label = label,
         linewidth = 2)

# ERROR BARS

if opts.mean:

    def one_sigma(points, bin_size, x, y):
        mean = []
        std = []
        for point in points:
            index = ((x >= point - bin_size) & \
                    (x <= point + bin_size))
            mean.append(np.mean(y[index]))
            std.append(np.std(y[index]))
        return mean, std

    mean, sigma = np.array(one_sigma(x, 0.01, z_spec, z_phot))

    plt.plot(x, mean, '--', color = 'k', label = 'Mean', linewidth = 2)
    plt.plot(x, mean + sigma, '-.', color = 'k', linewidth = 2,
            label = r'1-$\sigma$')
    plt.plot(x, mean - sigma, '-.', color = 'k', linewidth = 2)
    plt.plot(x, mean + 2 * sigma, ':', color = 'k', linewidth = 2,
            label = r'2-$\sigma$')
    plt.plot(x, mean - 2 * sigma, ':', color = 'k', linewidth = 2)
    plt.fill_between(x, mean - sigma, mean + sigma, facecolor='gray',
                    alpha=0.2)
    plt.fill_between(x, mean - 2 * sigma, mean + 2 * sigma,
                    facecolor='gray', alpha=0.1)

# FINISH PLOT
    
plt.xlabel(r'z$_{spec}$', {'fontsize' : 20})
if opts.diff:
    plt.ylabel(r'z$_{spec}$ - z$_{phot}$', {'fontsize' : 20})
else:
    plt.ylabel(r'z$_{phot}$', {'fontsize' : 20})
plt.legend(loc = 'upper left')
if not opts.output_file:
    plt_name = opts.input_file + '_analyz2.pdf'
else:
    plt_name = opts.output_file + '.pdf'
plt.savefig(plt_name)
print "Plot saved to:", plt_name

# PHOTOZ ERROR PLOT

if opts.zphoterr_col:

    plt.figure(2)
    hist2, xedges, yedges = np.histogram2d(z_phot, z_phot_err,
                                           bins = opts.n_bins)
    extent2 = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    hist2[hist2 == 0.0] = None
    plt.imshow(hist2.transpose()[::-1], extent = extent2,
           interpolation = opts.interp, cmap = cm.jet,
           aspect = 'auto', norm = LogNorm())
    plt.colorbar()
    plt.autoscale(False)
    if not opts.x_range:
        plt.xlim(min(z_spec), max(z_spec))
    else:
        plt.xlim(opts.x_range)
    if not opts.y2_range:
        plt.ylim(0.0, 0.3)
    else:
        plt.ylim(opts.y2_range)
    plt.xlabel(r'z$_{phot}$', {'fontsize' : 20})
    plt.ylabel(r'$\sigma_{z_{phot}}$', {'fontsize' : 20})
    if not opts.output_file:
        plt_name = opts.input_file + '_analyz2_err.pdf'
    else:
        plt_name = opts.output_file + '_err.pdf'
    plt.savefig(plt_name)
    print "Error Plot saved to:", plt_name
