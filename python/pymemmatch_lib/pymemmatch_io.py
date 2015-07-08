#######################
#   PYMEMMATCH_IO     #
#######################
# Samuel Farrens 2015 #
#######################

import numpy as np
from numpy.lib.recfunctions import append_fields

def read_mock(opts):

    """
    Function to read mock members file.
    """
    
    print 'Reading file:', opts.input_files[1],

    # read file        
    mock = np.genfromtxt(opts.input_files[1], dtype = "S", unpack = True)[opts.mock_cols]
    
    dtypes = [('m_mem_id', 'S22'), ('m_id', 'S22'), ('m_mass', 'f8'), ('m_ngal', 'int')]
        
    mock = np.core.records.fromarrays(mock, dtype = dtypes)

    # sort by member ID
    index = np.argsort(mock.m_mem_id)
        
    mock = mock[index]
    
    print '\tComplete:', len(mock)

    return mock

def read_obs(opts):

    """
    Function to read observed cluster members file.
    """

    print 'Reading file:', opts.input_files[0],

    # read file 
    obs = np.genfromtxt(opts.input_files[0], dtype="S", unpack = True)[opts.obs_cols]
    
    dtypes = [('id', 'i4'), ('proxy', 'int'), ('z', 'f8'), ('mem_id', 'S22')]
    
    obs = np.core.records.fromarrays(obs, dtype = dtypes)

    # sort by member ID
    index = np.argsort(obs.mem_id)

    obs = obs[index]
        
    print '\tComplete:', len(obs)

    return obs

def print_matrix(matrix, x_range, output_name):

    """
    Function to output mass-observable matrix.
    """

    file_name = output_name + '.pycy.matrix.txt'
    
    output = np.fliplr(np.transpose(np.vstack([matrix, x_range])))

    np.savetxt(file_name, output, fmt = '%.3f')
    
    print 'Data saved to:', file_name
