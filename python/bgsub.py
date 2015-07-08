#! /Users/Rowen/Documents/Library/anaconda/bin/python

## @file bgsub.py
#
#  BGSUB
#
#  Script for performing background subtraction
#  of H-alpha emitting galaxies.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
from bgsub_lib import *
from functions.interface import h_line

##
#  Code Main.
def main():

    # Get code arguments
    opts = bgsub_opts.get_opts()

    # Read input files
    h_line()
    
    c_data = bgsub_io.read_clusters(opts)
    g_data = bgsub_io.read_galaxies(opts)
    m_data = None
    
    if len(opts.input_file) > 2:
        m_data = bgsub_io.read_matches(opts)

    # Subtract H-alpha background
    h_line()

    if opts.c_id:
        bgsub_main.single(c_data, g_data, opts)

    elif opts.l_lim:
        bgsub_main.stack(c_data, g_data, opts, m_data)

    h_line()
    
if __name__ == "__main__":
    main()
