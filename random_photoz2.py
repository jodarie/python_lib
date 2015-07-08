#! /Users/Rowen/Documents/Library/anaconda/bin/python

#####################
# RANDOM PHOTOZ v.2 #
#####################

import argparse, os.path
import numpy as np
from scipy import interpolate

# GENERATE GAUSSIAN SPREAD

def gaussian(mu, sigma):
    return np.random.normal(mu, sigma)

# GENERATE OUTLIERS

def outliers(var, frac, z_min, z_max):
    index = np.arange(len(var))
    np.random.shuffle(index)
    index = index[:int(len(var) * frac / 100.)]
    var[index] = np.random.random(len(var[index])) * \
        (z_max - z_min) + z_min
    return var

# READ ARGUMENTS

parser = argparse.ArgumentParser("RANDOM PHOTOZ 2 OPTIONS")
parser.add_argument("-i", "--input_file",
                  dest = "input_file",
                  help = "Input file name.")
parser.add_argument("-c", "--col", dest = "col",
                  type = int,
                  help = "Redshift column number.")
parser.add_argument("-s", "--sigma_file", dest = "sigma_file",
                  help = "Standard deviation file.")
parser.add_argument("-f", "--fraction", dest = "frac",
                  type = float, default = 0.0,
                  help = "Outlier fraction.")
opts = parser.parse_args() 

if not opts.input_file:
    parser.err('Error: Input file not specified.')
 
if not os.path.exists(opts.input_file):
    parser.err('Error: Input file not found.')

if opts.col < 1:
    parser.err('Error: Column number must be > 0.')

if not opts.sigma_file:
    parser.err('Error: Sigma not specified.')

if not os.path.exists(opts.sigma_file):
    parser.err('Error: Sigma file not found.')

if opts.frac < 0.0 or opts.frac > 100.0:
    parser.err('Error: Fraction must be between 0.0 and 100.0.')

# READ FILES

data1 = np.genfromtxt(opts.input_file, unpack = True,
                     dtype = 'S')
z = np.array(data1[opts.col - 1, :], dtype = 'float')

sigma = np.genfromtxt(opts.sigma_file, unpack = True,
                      dtype = 'float')

# INTERPOLATE SIGMA VALUES

sig_func = interpolate.interp1d(sigma[0], sigma[2])

# GENERATE GAUSSIAN PHOTOMETRIC REDSHIFT

z_phot = np.random.normal(z, sig_func(z))

z_phot = outliers(z_phot, opts.frac, min(z), max(z))

# INSERT NEW COLUMN INTO DATA

new_data = np.insert(data1.transpose(), opts.col,
                     np.array([z_phot, sig_func(z)], dtype = 'S'),
                     axis = 1)

# SAVE NEW ARRAY TO FILE

out_file = opts.input_file + '_rand_photoz.dat'
np.savetxt(out_file, new_data, fmt = '%s')
print 'Saving array contents to:', out_file
