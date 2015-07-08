## @file pycymatch_match.py
#
#  PYCYMATCH MATCH METHODS 
#
#  Functions for finding matches
#  between observed and mock
#  clusters.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import sys
sys.path.append("/Users/Rowen/Documents/Codes/Python")

import numpy as np
import itertools
from functions import astro, comp
from pycymatch_io import print_matches

##
#  Find observed clusters that satisfy the
#  basic match criteria.
#
#  @param[in] i: Index of observed cluster.
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] dz: Redshift matching treshold.
#
#  @return List of indices of mock haloes
#  that satisfy the matching criteria for a
#  given observed cluster.
#
def basic_match(i, mock, obs, dz):

    return np.where((np.fabs(mock.z - obs.z[i]) <= \
                     (dz * (1 + mock.z))) &
                     (obs.ra[i] >= mock.minra) &
                     (obs.ra[i] <= mock.maxra) &
                     (obs.dec[i] >= mock.mindec) &
                     (obs.dec[i] <= mock.maxdec))[0]

##
#  Find observed clusters within r200 of mock
#  halo centre.
#
#  @param[in] i: Index of observed cluster.
#  @param[in] index: List of indices of matches.
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#
#  @return Indices of mathces and corresponding
#  probabilities of mock haloes to a given
#  observed cluster.
#
def r200_match(i, index, mock, obs):

    dists = astro.ang_sep((mock.ra[index], mock.dec[index]), (obs.ra[i], obs.dec[i])) * 60.0

    new_index = np.where((dists - (1.0 * mock.r200[index])) <= 0.0)[0]

    if len(new_index) > 0:
        dist_data = 1.0 - (dists[new_index] / mock.r200[index[new_index]]) ** 2  
        return index[new_index], comp.scale(dist_data, 0.0, np.sum(dist_data))
    else:
        return [], []

##
#  Find cylindrical matches between mock
#  haloes and observed clusters.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] dz: Redshift matching treshold.
#
#  @return Indices of mathces and corresponding
#  probabilities of mock haloes to observed
#  clusters.
#
def find_matches(mock, obs, dz):

    index_list = []
    dist_list = []
    
    for i in range(obs.size):
        index = basic_match(i, mock, obs, dz)
                        
        if len(index) > 0:
            x, y = r200_match(i, index, mock, obs)                              
            index_list.append(x)
            dist_list.append(y)
            
        else:
            index_list.append([])
            dist_list.append([])
            
    return index_list, dist_list

##
#  Produce mass-observable matrix.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] opts: List of arguments.
#
#  @return Mass-observable matrix, histogram
#  of all mock halo masses, histogram of
#  matched halo masses, x-range values of
#  halo masses, x-range values of observed
#  cluster mass proxies.
#
def mo_matrix(mock, obs, opts):

    n_mass_bins = comp.num_bins(opts.mass_bin[0], opts.mass_bin[1], opts.mass_bin[2])
    n_proxy_bins = comp.num_bins(opts.proxy_bin[0], opts.proxy_bin[1], opts.proxy_bin[2])
    
    mass_index = np.floor((mock.mass - opts.mass_bin[0]) / opts.mass_bin[2]).astype('int')
    proxy_index = np.floor((obs.proxy - opts.proxy_bin[0]) / opts.proxy_bin[2]).astype('int')

    mass_hist = np.zeros(n_mass_bins)
    for i in mass_index:
        mass_hist[i] += 1
    
    mass_x = comp.x_vals(n_mass_bins, opts.mass_bin[0], opts.mass_bin[2])
    proxy_x = comp.x_vals(n_proxy_bins, opts.proxy_bin[0], opts.proxy_bin[2])

    matrix = np.zeros((n_proxy_bins, n_mass_bins))
    
    index, dists = find_matches(mock, obs, 2.0 * opts.delta_z)

    tag_index = halo_tag(mock, np.copy(index))

    mass_hist2 = np.zeros(n_mass_bins)
    for i in mass_index[tag_index]:
        mass_hist2[i] += 1

    print_matches(mock, obs, index, opts)

    for i in range(obs.size):
        for j in range(len(index[i])):
            matrix[proxy_index[i], mass_index[index[i][j]]] += dists[i][j]
            
    return np.nan_to_num(np.flipud(matrix)), mass_hist, mass_hist2, mass_x, proxy_x

##
#  Function to scale matrix lines to 1.0.
#
#  @param[in] matrix: Input matrix.
#
#  @return Scaled matrix.
#
def scale_matrix(matrix):
    
    for i in range(len(matrix)):
        matrix[i] = comp.scale(matrix[i], 0.0, np.sum(matrix[i]))

    return matrix

##
#  Function to tag matched Haloes.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] index: List of indices of matches.
#
#  @return List of indeces of unique
#  matches to mock haloes.
#
def halo_tag(mock, index):

    new_index = np.unique(list(itertools.chain.from_iterable(index.tolist())))

    print 'Found', len(new_index), 'out of', len(mock), 'haloes.', \
      '%.2f%%' % ((float(len(new_index)) / float(len(mock))) * 100.0)

    return new_index
