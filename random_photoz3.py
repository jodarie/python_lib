#! /Users/Rowen/Documents/Library/anaconda/bin/python

#####################
# RANDOM PHOTOZ v.3 #
#####################

import argparse, os.path
import numpy as np
from functions import comp
import matplotlib.pyplot as plt

# READ ARGUMENTS

parser = argparse.ArgumentParser("RANDOM PHOTOZ 3 OPTIONS")
parser.add_argument("-i", "--input_file",
                  dest = "input_file",
                  help = "Input file name.")
parser.add_argument("-c", "--col", dest = "col",
                  type = int,
                  help = "Redshift column number.")
parser.add_argument("-p", "--pdf_file", dest = "pdf_file",
                  help = "PDF deviation file.")
opts = parser.parse_args()

if not opts.input_file:
    parser.err('Error: Input file not specified.')
 
if not os.path.exists(opts.input_file):
    parser.err('Error: Input file not found.')

if opts.col < 1:
    parser.err('Error: Column number must be > 0.')

if not opts.pdf_file:
    parser.err('Error: PDF not specified.')

if not os.path.exists(opts.pdf_file):
    parser.err('Error: PDF file not found.')

# READ FILES

data1 = np.genfromtxt(opts.input_file, unpack = True,
                     dtype = 'S')
z_spec = np.array(data1[opts.col - 1, :], dtype = 'float')

pdf = np.genfromtxt(opts.pdf_file, unpack = True,
                      dtype = 'float')

n_bins = pdf.shape[1]

bin_size = np.abs(pdf[0, 1] - pdf[0, 0])
z_min = np.min(pdf[0])
z_max = np.max(pdf[0]) + bin_size

index = np.array(np.floor(np.round((z_spec - z_min) / bin_size, 8)), dtype = 'int')

# GENERATE Z_PHOTS THAT FOLLOW THE PDF

z_rand = []
z_err = []
for i in range(len(z_spec)):
    z_rand.extend(np.random.choice(pdf[0], 1, p = pdf[index[i] + 1]))
    z_err.extend(np.random.choice(pdf[0], 1, p = pdf[index[i] + n_bins + 1]))
    z_err[i] = np.abs(z_err[i] + np.random.uniform(-0.01 / 2.0, 0.01 / 2.0))
    z_rand[i] = np.abs(z_rand[i] + np.random.uniform(-z_err[i] / 2.0, z_err[i] / 2.0))
z_rand = np.array(z_rand)
z_err = np.array(z_err)

# INSERT NEW COLUMN INTO DATA

new_data = np.insert(data1.transpose(), opts.col,
                     np.array([z_rand, z_err], dtype = 'S'),
                     axis = 1)

# SAVE NEW ARRAY TO FILE

out_file = opts.input_file + '_rand_photoz.dat'
np.savetxt(out_file, new_data, fmt = '%s')
print 'Saving array contents to:', out_file
