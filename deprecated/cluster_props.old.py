#! /Users/Rowen/Documents/Library/anaconda/bin/python

######################
# CLUSTER_PROPS V1.0 #
######################

import sys, argparse, os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from functions import cosmo

# READ ARGUMENTS

parser = argparse.ArgumentParser('CLUSTER_PROPS OPTIONS:')
parser.add_argument('-i', '--input_file', dest = 'input_file',
                    help = 'Input file name.')
parser.add_argument('-o', '--output_file', dest = 'output_file',
                    help = 'Output file name.')
parser.add_argument('-c', '--columns', dest = 'columns',
                    default = np.arange(1, 12),
                    nargs = '+', type = int, help = 'Column numbers of' + \
                    'cluster properties: ID RA RA_err Dec Dec_err z z_err' + \
                    'N_gal S/N R(arcmin) A(arcmin^2) ' + \
                    '[Default = 1 2 3 4 5 6 7 8 9 10 11]')
opts = parser.parse_args()

if not opts.input_file:
    print 'Error: Input file not specified.'
    exit()
 
if not os.path.exists(opts.input_file):
    print 'Error: Input file not found.'
    exit()

# READ FILE

data = np.genfromtxt(opts.input_file, unpack = True,
                     dtype = 'float')

r_phys = (data[9] * (np.pi/(180.0 * 60.0))) * (3000) * \
  map(cosmo.angdidis, data[5], np.ones(len(data[5])) * \
      0.3, np.ones(len(data[5])) * 0.7)

a_phys = np.pi * r_phys

# MAKE PLOTS

def histo(name, x, bins, xran, yran, xlog, ylog, xlab, ylab, tit):
    hist, bins = np.histogram(x, bins = bins)
    plt.figure()
    plt.plot(bins, hist, 'b-')
    if xran:
        plt.xlim(xran)
    if yran:
        plt.ylim(yran)
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(tit)
    plt.savefig(name)
    
def d_map(name, x, y, bins, interp, norm, xran, yran, xlog, ylog, xlab, ylab, tit):
    hist, xedges, yedges = np.histogram2d(x, y, bins = bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.figure()
    plt.imshow(hist.transpose()[::-1], extent = extent,
               interpolation = interp, cmap = cm.jet,
               aspect = 'auto', norm = norm)
    plt.colorbar()
    plt.autoscale(False)
    if xran:
        plt.xlim(xran)
    if yran:
        plt.ylim(yran)
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(tit)
    plt.savefig(name)

d_map('cluster_dist.pdf', data[1], data[3], 30, 'spline16', None,
      False, False, False, False, 'RA', 'Dec', 'Cluster Distribution')

d_map('cluster_ngal.pdf', data[5], data[7], 100, 'spline16', LogNorm(),
      False, False, False, False, 'z', r'N$_{gal}$', 'Cluster Members')

d_map('cluster_logngal.pdf', data[5], np.log10(data[7]), 50, 'spline16',
      LogNorm(), False, False, False, False, 'z', r'$\log_{10}$N$_{gal}$',
      'Cluster Members')

d_map('cluster_sn.pdf', data[5], data[8], 100, 'spline16', LogNorm(),
      False, False, False, False, 'z', 'S/N',
      'Cluster Signa-to-Noise Ratios')

d_map('cluster_logsn.pdf', data[5], np.log10(data[8]), 50, 'spline16',
      LogNorm(), False, False, False, False, 'z', r'$\log_{10}$S/N',
      'Cluster Signa-to-Noise Ratios')

d_map('cluster_size_arcmin.pdf', data[5], data[9], 100, 'spline16', LogNorm(),
      False, False, False, False, 'z', 'R [arcmin]', 'Cluster Sizes')

d_map('cluster_area_armin2.pdf', data[5], data[10], 150, 'spline16', LogNorm(),
      False, False, False, False, 'z', r'A [arcmin$^2$]',
      'Cluster Areas')

d_map('cluster_size_mpc.pdf', data[5], r_phys, 100, 'spline16', LogNorm(),
      False, False, False, False, 'z', 'R [Mpc h$^{-1}$]', 'Cluster Sizes')

d_map('cluster_area_mpc2.pdf', data[5], a_phys, 100, 'spline16', LogNorm(),
      False, False, False, False, 'z', 'R [Mpc$^2$ h$^{-2}$]', 'Cluster Areas')

