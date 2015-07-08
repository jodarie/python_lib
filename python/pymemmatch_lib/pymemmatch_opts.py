#######################
#   PYMEMMATCH_OPTS   #
#######################
# Samuel Farrens 2015 #
#######################

import argparse

def get_opts():

    """
    Function to read in code arguments.
    """

    m_col_help = ('Set column numbers for mock properties: MEM_ID, HALO_ID, HALO_MASS, HALO_NGAL.' + \
                  '[Default: 1 6 7 9]')
                
    o_col_help = ('Set column numbers for cluster properties: OBS_ID, OBS_PROXY, OBS_Z, MEM_ID. ' + \
                  '[Default 1 2 3 4]')

    parser = argparse.ArgumentParser('PYMEMMATCH OPTIONS:')
    
    parser.add_argument('-i', '--input_files', dest = 'input_files', nargs = '+',
                        help = 'Input file names. (1-Observerd Cluster Members, 2-Mock Halo Members)')
    
    parser.add_argument('--proxy_bin', dest = 'proxy_bin', default = [0.0, 3.0, 0.1],
                        nargs = '+', type = float, help = 'Mass proxy bin values.' + \
                        ' [Default: 0.0, 3.0, 0.2]')
                        
    parser.add_argument('--mass_bin', dest = 'mass_bin', default = [13.0, 15.5, 0.1],
                        nargs = '+', type = float, help = 'Mass bin values.' + \
                        ' [Default: 13.0, 15.5, 0.1]')
                        
    parser.add_argument('--z_bin', dest = 'z_bin', default = [0.0, 3.0, 0.1],
                        nargs = '+', type = float, help = 'Redshift bin values.' + \
                        ' [Default: 0.0, 3.0, 0.1]')
    
    parser.add_argument('--mock_cols', dest = 'mock_cols', default = [0, 5, 6, 8],
                        nargs = '+', type = int, help = m_col_help )
    
    parser.add_argument('--obs_cols', dest = 'obs_cols', default = range(0, 4),
                        nargs = '+', type = int, help = o_col_help )
    
    check_opts(parser)
    
    return parser
    

def check_opts(parser):

    """
    Function to check parser arguments.
    """

    opts = parser.parse_args()
    
    if not opts.input_files:
        parser.error('argument --input_files: file names not provided')
    
    if not len(opts.input_files) == 2:
        parser.error('argument --input_files: requires 2 input file names\n' \
                    + 'e.g. --input_files obs_cluster_file mock_halo_file')
                 
    if not len(opts.proxy_bin) == 3:
        parser.error('argument --proxy_bin: requires 3 input values\n' \
                    + 'e.g. --proxy_bin min_value max_value bin_size')
                 
    if not len(opts.mass_bin) == 3:
        parser.error('argument --mass_bin: requires 3 input values\n' \
                    + 'e.g. --mass_bin min_value max_value bin_size')
                 
    if not len(opts.z_bin) == 3:
        parser.error('argument --z_bin: requires 3 input values\n' \
                    + 'e.g. --z_bin min_value max_value bin_size')
                           
    if not len(opts.mock_cols) == 4:
        parser.error('argument --mock_cols: requires 4 input values')
            
    if not len(opts.obs_cols) == 4:
        parser.error('argument --obs_cols: requires 4 input values')    
