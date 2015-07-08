#! /Users/Rowen/Documents/Library/anaconda/bin/python

## @file pycymatch_v5.py
#
#  PYCYMATCH 
#
#  Script for performing 
#  cylindrical matching between
#  observed and mock clusters.
#  
#  @author Samuel Farrens
#  @version 5.1
#  @date 2015
#

import numpy as np

from functions import errors
from pycymatch_lib.pycymatch_opts import get_opts
from pycymatch_lib.pycymatch_extra import h_line
from pycymatch_lib.pycymatch_match2 import *
from pycymatch_lib.pycymatch_cp import *
from pycymatch_lib.pycymatch_fit import get_fits
from pycymatch_lib.pycymatch_io import *
from pycymatch_lib.pycymatch_plot import *

##
#  Code Main.
def main():

    ##
    # Read arguments
    opts = get_opts().parse_args()

    ##
    # Check input files
    for file in opts.input_files:
        errors.file_name_error(file)

    h_line()
    
    ##
    # Read mock halo catalogue
    mock = read_mock(opts)

    ##
    # Read observed catalogue
    obs = read_obs(opts)

    h_line()

    ##
    # Find Matches
    matches = find_matches(mock, obs, 2.0 * opts.delta_z, opts)
    print_matches(mock[matches[2]], obs, matches[0], matches[1], opts)

    ##
    # Define completeness and purity of sample
    c_matrix = get_completeness(mock, matches, opts)
    p_matrix = get_purity(obs, matches, opts)

    ##
    # Define mass-observable matrix for matched objects
    matrix, hm_matrix, ranges = mo_matrix(mock, obs, matches, opts)
        
    h_line()
    
    ##
    # Make plots
    make_plots(matrix, hm_matrix, ranges, opts)
    plot_complete(c_matrix[0], 'mass', opts)
    plot_complete(c_matrix[1], 'ngal', opts)
    plot_pure(p_matrix, opts)
    
    ##
    # Save matrix to file
    print_matrix(matrix, hm_matrix, ranges, opts)

    h_line()

if __name__ == "__main__":
    main()
