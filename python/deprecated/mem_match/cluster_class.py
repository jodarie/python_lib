#################
# CLUSTER CLASS #
#################

import numpy as np
import functions

class Clusterx:
    
    '''
    Class for storing cluster properties.
    '''

    def __init__(self, c_id):
        '''
        Function that initialises a Cluster instance.
        '''
        self.id = c_id
        self.mem = []
        self.n_mem = 0

    def add_mem(self, member):
        '''
        Function to add a galaxy member to a Cluster instance.
        '''
        self.mem.append(member)
        self.count_mem()

    def count_mem(self):
        '''
        Function to count number of members.
        '''
        self.n_mem = len(self.mem)

    def props(self, proxy_col, z_col):
        '''
        Function to determine cluster properties.
        '''
        self.proxy = self.mem[0][proxy_col - 1]
        self.z = self.mem[0][z_col - 1]
        
    def halo_count(self):
        '''
        Function to count the number of mock haloes with
        common members.
        '''
        mem_mock_ids = []
        self.mem_mass = []
        for member in self.mem:
            mem_mock_ids.append(member[-2])
            self.mem_mass.append(member[-3])
        self.mock_id = np.unique(mem_mock_ids, return_counts = True)
        self.mass = np.unique(self.mem_mass, return_counts = True, return_index = True)
    
    def mass_hist(self, min_val, max_val, bin_size):
        '''
        Function to produce mass histogram.
        '''
        self.hist_x, self.hist = functions.get_hist(self.mem_mass,
                                                    min_val, max_val, bin_size)
        
