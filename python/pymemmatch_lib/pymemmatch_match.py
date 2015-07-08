#######################
#   PYMEMMATCH_MATCH  #
#######################
# Samuel Farrens 2015 #
#######################

import numpy as np
from functions import comp
from pymemmatch_cc import Clusterx

def find_matches(mock, obs, opts):

    """
    Function to find matching galaxy members between mock haloes
    and observed clusters.
    """
    
    obs = obs[np.in1d(obs.mem_id, mock.m_mem_id, assume_unique = True)]
    
    mock = mock[np.in1d(mock.m_mem_id, obs.mem_id, assume_unique = True)]
    
    merged = np.lib.recfunctions.merge_arrays([obs, mock], flatten = True,
                                              usemask = False)

    clusters = []
    count = 0

    for id_val in np.unique(obs.id):
        clusters.append(Clusterx(count))
         
        for member in merged[obs.id == id_val]:
            clusters[count].add_mem(member)
            
        count += 1
    
    for cluster in clusters:
        cluster.props()
        cluster.halo_count()
        cluster.mass_hist(opts.mass_bin)

    return clusters

def mo_matrix(mock, obs, opts):

    """
    Function to determine the mass-observable matrix.
    """

    n_mass_bins = comp.num_bins(opts.mass_bin[0], opts.mass_bin[1], opts.mass_bin[2])
    n_proxy_bins = comp.num_bins(opts.proxy_bin[0], opts.proxy_bin[1], opts.proxy_bin[2])

    mass_x = comp.x_vals(n_mass_bins, opts.mass_bin[0], opts.mass_bin[2])
    proxy_x = comp.x_vals(n_proxy_bins, opts.proxy_bin[0], opts.proxy_bin[2])
    
    matrix = np.zeros((n_proxy_bins, n_mass_bins))

    clusters = find_matches(mock, obs, opts)

    for cluster in clusters:
    
        proxy_bin = n_proxy_bins - 1 - \
          comp.find_bin(np.log10(cluster.proxy), opts.proxy_bin[0], opts.proxy_bin[2])

        for i in range(n_proxy_bins):
            if proxy_bin == i and opts.z_bin[0] <= cluster.z < opts.z_bin[1]:
                matrix[proxy_bin] += comp.scale(cluster.hist, min(cluster.hist),
                                                np.sum(cluster.hist))
                   
    for i in range(n_proxy_bins):
        matrix[i] = comp.scale(matrix[i], min(matrix[i]), np.sum(matrix[i]))

    return matrix, mass_x, proxy_x
