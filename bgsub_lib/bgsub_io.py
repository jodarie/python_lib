## @file bgsub_io.py
#
#  BGSUB IO FUNCTIONS 
#
#  Functions for managing
#  input and output
#  operations.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
import pandas as pd
import random

##
#  Read ascii file of H-alpha galaxies
#  and return pandas DataFrame.
#
#  @param[in] opts: List of arguments.
#
#  @return DataFrame of galaxies.
#
def read_galaxies(opts):

    data = np.genfromtxt(opts.input_file[0], unpack = True,
                         dtype = 'S', usecols = (np.array(opts.g_cols) - 1))

    print ' Successfully read:', opts.input_file[0], '(%i galaxies)' % data.shape[1]

    index = random.sample(range(len(data[0])), int(float(len(data[0])) * (opts.complete / 100.0)))
    print ' - Using {0:.1f}% ({1:d} galaxies)'.format(opts.complete, len(index))
        
    return pd.DataFrame({'id' : data[0, index],
                         'ra' : np.array(data[1, index], dtype = 'float64'),
                         'dec' : np.array(data[2, index], dtype = 'float64'),
                         'z' : np.array(data[3, index], dtype = 'float64')})

##
#  Read ascii file of cluster candidates
#  and return pandas DataFrame.
#
#  @param[in] opts: List of arguments.
#
#  @return DataFrame of clusters.
#
def read_clusters(opts):

    data = np.genfromtxt(opts.input_file[1], unpack = True,
                         dtype = 'S', usecols = (np.array(opts.c_cols) - 1))

    print ' Successfully read:', opts.input_file[1], '(%i clusters)' % data.shape[1]

    return pd.DataFrame({'id' : data[0],
                         'ra' : np.array(data[1], dtype = 'float64'),
                         'dec' : np.array(data[2], dtype = 'float64'),
                         'z' : np.array(data[3], dtype = 'float64'),
                         'lamb' : np.log10(np.array(data[4], dtype = 'float64'))})

##
#  Read ascii file of cluster matches
#  to mock haloes.
#
#  @param[in] opts: List of arguments.
#
#  @return Dictionary of matches.
#
def read_matches(opts):

    data = np.genfromtxt(opts.input_file[2], dtype = 'S',
                         usecols = (0, 6))

    return dict(data)
