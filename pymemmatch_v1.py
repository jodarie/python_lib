#! /Users/Rowen/Documents/Library/anaconda/bin/python

#######################
#   PYMEMMATCH V.5.0  #
#######################
# Samuel Farrens 2015 #
#######################

import numpy as np

from functions import errors
from pymemmatch_lib.pymemmatch_opts import get_opts
from pycymatch_lib.pycymatch_extra import h_line
from pymemmatch_lib.pymemmatch_match import mo_matrix
from pymemmatch_lib.pymemmatch_io import *
from pycymatch_lib.pycymatch_plot import *

# READ ARGUMENTS 

opts = get_opts().parse_args()

# READ DATA 

# - check files
for file in opts.input_files:
    errors.file_name_error(file)

h_line()

# - read mock catalogue
mock = read_mock(opts)

# - read observed catalogue
obs = read_obs(opts)

h_line()

# FIND MATCHES

matrix, x_range, z_range = mo_matrix(mock, obs, opts)

# PLOT

mass_obs(matrix, opts)
p_lambda(matrix, x_range, z_range, opts)

# SAVE MATRIX TO FILE

print_matrix(matrix, x_range, opts.input_files[0])

h_line()
