#! /Users/Rowen/Documents/Library/anaconda/bin/python

#################
# RANDOM PHOTOZ #
#################

import sys, optparse, os.path
import numpy as np
import matplotlib.pyplot as plt

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

parser = optparse.OptionParser()
parser.add_option("-i", "--input_file",
                  dest = "input_file",
                  help = "Input file name.")
parser.add_option("-c", "--col", dest = "col",
                  type = "int",
                  help = "Redshift column number.")
parser.add_option("-s", "--sigma", dest = "sigma",
                  type = "float",
                  help = "Standard deviation.")
parser.add_option("-f", "--fraction", dest = "frac",
                  type = "float",
                  help = "Outlier fraction.")
(opts, args) = parser.parse_args()

if not opts.input_file:
    print 'Error: Input file not specified.'
    exit()
 
if not os.path.exists(opts.input_file):
    print 'Error: Input file not found.'
    exit()

if opts.col < 1:
    print 'Error: Column number must be > 0.'
    exit()

if opts.sigma < 0.0:
    print 'Error: Sigma must be >= 0.0.'
    exit()

if opts.frac < 0.0 or opts.frac > 100.0:
    print 'Error: Fraction must be between 0.0 and 100.0.'
    exit()

# READ FILE

data = np.genfromtxt(opts.input_file, unpack = True,
                     dtype = 'S')
z = np.array(data[opts.col - 1, :], dtype = 'float')

# GENERATE GAUSSIAN PHOTOMETRIC REDSHIFT

z_phot = outliers(gaussian(z, opts.sigma * (1 + z)), opts.frac,
                  min(z), max(z))

# PLOT Z VS Z_PHOT

print 'Close plot to print results to file.'

plt.hist2d(z, z_phot, 200)
plt.xlabel(r'z')
plt.ylabel(r'z$_{phot}$')
plt.title('Gaussian Photometric Redshifts with ' + \
          str(opts.frac) + '% outliers.')
plt.xlim(min(z), max(z))
plt.ylim(min(z), max(z))
plt.show()

# INSERT NEW COLUMN INTO DATA

new_data = np.insert(data.transpose(), opts.col,
                     np.array(z_phot, dtype = 'S'),
                     axis = 1)

# SAVE NEW ARRAY TO FILE

out_file = opts.input_file + '_rand_photoz.dat'
np.savetxt(out_file, new_data, fmt = '%s')
print 'Saving array contents to:', out_file
