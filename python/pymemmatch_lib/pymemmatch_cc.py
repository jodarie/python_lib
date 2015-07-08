#######################
#   PYMEMMATCH_CC     #
#######################
# Samuel Farrens 2015 #
#######################

import sys
sys.path.append("/Users/Rowen/Documents/Codes/Python")

import numpy as np
from functions import comp

def get_hist(data, bin_vals):
    
    """
    Generate histogram.
    """
    
    n_bins = comp.num_bins(bin_vals[0], bin_vals[1], bin_vals[2])
    hist = np.histogram(data, bins = n_bins,
                        range = (bin_vals[0], bin_vals[1]))
    return np.array(hist[0], dtype = 'float')

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

    def props(self):
        
        '''
        Function to determine cluster properties.
        '''
        
        self.proxy = self.mem[0][1]
        self.z = self.mem[0][2]
        
    def halo_count(self):
        
        '''
        Function to count the number of mock haloes with
        common members.
        '''
        
        mem_mock_ids = []
        self.mem_mass = []
        for member in self.mem:
            mem_mock_ids.append(member[-3])
            self.mem_mass.append(member[-2])
        self.mock_id = np.unique(mem_mock_ids, return_counts = True)
        self.mass = np.unique(self.mem_mass, return_counts = True,
                              return_index = True)
    
    def mass_hist(self, bin_vals):
        
        '''
        Function to produce mass histogram.
        '''
        
        self.hist = get_hist(self.mem_mass, bin_vals)
        
