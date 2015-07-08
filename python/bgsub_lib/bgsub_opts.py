## @file bgsub_opts.py
#
#  BGSUB OPTION FUNCTIONS 
#
#  Functions for retreiving
#  code arguments.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import argparse

##
#  Function to read in code arguments.
#
#  @return List of arguments.
#
def get_opts():

    c_col_help = ('Set column numbers for cluster properties: ID, RA, DEC, Z, LAMBDA' + \
                ' [Default 1 2 3 4 5]')

    g_col_help = ('Set column numbers for galaxy properties: ID, RA, DEC, Z' + \
                ' [Default 1 2 3 4]')

    parser = argparse.ArgumentParser('BGSUB OPTIONS:',
                                     formatter_class = argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input_file', dest = 'input_file', nargs = '+',
                        help = 'Input file name.')

    parser.add_argument('--cluster_cols', dest = 'c_cols', default = range(1, 6),
                        nargs = '+', type = int, help = c_col_help)
        
    parser.add_argument('--galaxy_cols', dest = 'g_cols', default = range(1, 5),
                        nargs = '+', type = int, help = g_col_help)

    parser.add_argument('-n', '--cluster_id', dest = 'c_id', default = '1',
                        help = 'Cluster ID. [Default: 1]')

    parser.add_argument('--ra_lim', dest = 'ra_lim', default = (45.0, 55.0),
                        nargs = '+', type = float, help = 'RA Limits [Default: 45.0 55.0]')

    parser.add_argument('--dec_lim', dest = 'dec_lim', default = (-5.0, 5.0),
                        nargs = '+', type = float, help = 'Dec Limits [Default: -5.0 5.0]')

    parser.add_argument('--z_lim', dest = 'z_lim', default = (0.9, 1.8),
                        nargs = '+', type = float, help = 'Redshift Limits [Default: 0.9 1.8]')
    
    parser.add_argument('--l_lim', dest = 'l_lim', nargs = '+', type = float, help = 'Lambda Limits')

    parser.add_argument('--n_bins', dest = 'n_bins', default = 18, type = int,
                        help = 'Number of redshift bins. [Default: 18]')

    parser.add_argument('--n_rand', dest = 'n_rand', default = 1000, type = int,
                        help = 'Number of random background samples. [Default: 1000]')

    parser.add_argument('-c', '--completeness', dest = 'complete', default = 100.0, type = float,
                        help = 'Completeness of H-alpha catalogue. [Default: 100.0]')
    
    opts = check_opts(parser)
    
    return opts

##
#  Function to read in code arguments.
#
#  @param[in] parser: List of arguments.
#
def check_opts(parser):

    opts = parser.parse_args()
    
    if not opts.input_file:
        parser.error('argument --input_file: file name not provided.')

    if opts.l_lim:
        opts.c_id = None

    if opts.complete <= 0.0 or opts.complete > 100.0:
        raise TypeError('Invalid value [{0:.1f}].'.format(opts.complete) + \
                        'Completeness must be in range 0.0 < C <= 100.0.')

    return opts
