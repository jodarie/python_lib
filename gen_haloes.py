#! /Users/Rowen/Documents/Library/anaconda/bin/python

## @file gen_haloes.py
#
#  GENERATE HALOES
#
#  Format halo catalogue for
#  matching codes.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import math, argparse
import numpy as np
from library import const
from functions import astro, cosmo2
from halo_methods.virial import r200

# READ ARGUMENTS

parser = argparse.ArgumentParser('GEN_HALOES OPTIONS:')
parser.add_argument('-i', '--input_file', dest = 'input_file',
                    help = 'Input file name.')
parser.add_argument('-c', '--columns', dest = 'columns', type = int,
                    nargs = '+', help = 'File column numbers: Halo_ID Mass Central RA Dec Z')
parser.add_argument('-m', '--mass_limit', dest = 'mass_limit',
                    default = 13.0, type = float, help = 'Mass limit. [Default: 13.0]')
parser.add_argument('--mem_limit', dest = 'mem_limit',
                    default = 5, type = int, help = 'Member limit. [Default: 5]')
opts = parser.parse_args()

if not opts.columns:
    raise ValueError('Column numbers must be specified.')

# READ DATA

data = np.genfromtxt(opts.input_file, dtype="S", unpack = True)

# DEFINE HALO ID DICTIONARY

index1 = np.array(data[opts.columns[1] - 1], dtype = 'float') >= opts.mass_limit

u_id = np.unique(data[opts.columns[0] - 1, index1], return_counts = True)

u_dict = dict(zip(*u_id))

# DEFINE CONSTANTS

H0 = 100 #km/s/Mpc

# PRINT HALO FILE INFO

for i in range(len(u_id[0])):
    
    x = (u_id[0][i] == data[opts.columns[0] - 1])

    id = data[opts.columns[0] - 1, x][0]
    central = data[opts.columns[2] - 1, x][0]
    ra = np.median(np.array(data[opts.columns[3] - 1, x], dtype = 'f'))
    min_ra = np.min(np.array(data[opts.columns[3] - 1, x], dtype = 'f'))
    max_ra = np.max(np.array(data[opts.columns[3] - 1, x], dtype = 'f'))
    dec = np.median(np.array(data[opts.columns[4] - 1, x], dtype = 'f'))
    min_dec = np.min(np.array(data[opts.columns[4] - 1, x], dtype = 'f'))
    max_dec = np.max(np.array(data[opts.columns[4] - 1, x], dtype = 'f'))
    z = np.median(np.array(data[opts.columns[5] - 1, x], dtype = 'f'))
    rich = u_dict[id]
    mass = np.array(data[opts.columns[1] - 1, x][0], dtype = 'float')
    r200_val = r200(mass, H0, 0.0, *const.BASE) #Mpc
    da = cosmo2.d_angdi(z, *const.BASE) * cosmo2.d_H(H0) #Mpc
    r200_val_arcmin = astro.rad2deg(r200_val / da) * 60.0 #arcmin
    
    if rich >= opts.mem_limit:
        print id, central, ra, dec, z, rich, mass, r200_val_arcmin, min_ra, max_ra, min_dec, max_dec
