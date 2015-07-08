## @file bgsub_bgsub.py
#
#  BGSUB FUNCTIONS 
#
#  Functions for subtracting
#  H-alpha galaxy background.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
from functions.astro import *
from functions.cosmo2 import *
from library import const

##
#  Calcualte distances between points
#  in Mpc/h.
#
#  @param[in] pos1: List of galaxy positions.
#  @param[in] pos2: Cluster position.
#  @param[in] z: Cluster redshift.
#  @param[in] H0: Hubble constant.
#
#  @return Arrays of distances in Mpc/h.
#
def get_dists(pos1, pos2, z, H0):

    return deg2rad(ang_sep((pos1[0], pos1[1]), (pos2[0], pos2[1]))) * \
      d_angdi(z, *const.BASE) * d_H(H0)

#
#  Calcualte H-alpha background level.
#
#  @param[in] gals: List of galaxies.
#  @param[in] opts: List of arguments.
#
#  @return Histogram of background level.
#
def get_bg(gals, opts):

    if opts.c_id:
        hist_sum = np.zeros(opts.n_bins)
    else:
        hist_sum = np.zeros(2 * opts.n_bins)

    for i in range(opts.n_rand):
        
        ra = np.random.uniform(*opts.ra_lim)
        dec = np.random.uniform(*opts.dec_lim)
        z = np.random.uniform(*opts.z_lim)

        dists = get_dists((gals.ra, gals.dec), (ra, dec), z, 70.0)
    
        index = np.where(dists <= 1.0)[0]

        if opts.c_id:
            data = gals.z[index]
            limits = opts.z_lim
            n_bins = opts.n_bins
            
        else:
            data = gals.z[index] - z
            x = opts.z_lim[1] - opts.z_lim[0]
            limits = (-x, x)
            n_bins = 2 * opts.n_bins
        
        hist_sum += np.histogram(data, range = limits, bins = n_bins)[0]

    return hist_sum / float(opts.n_rand)

#
#  Calcualte H-alpha background +
#  foreground level.
#
#  @param[in] gals: List of galaxies.
#  @param[in] clt: DataFrame element of
#  given cluster. 
#  @param[in] opts: List of arguments.
#
#  @return Histogram of background +
#  foreground level.
#
def get_bgfg(gals, clt, opts):

    dists = get_dists((gals.ra, gals.dec), (float(clt.ra), float(clt.dec)), float(clt.z), 70.0)
      
    index = np.where(dists <= 1.0)[0]

    if opts.c_id:
        data = gals.z[index]
        limits = opts.z_lim
        n_bins = opts.n_bins

    else:
        data = gals.z[index] - float(clt.z)
        x = opts.z_lim[1] - opts.z_lim[0]
        limits = (-x, x)
        n_bins = 2 * opts.n_bins

    hist = np.histogram(data, range = limits, bins = n_bins)
    
    return hist[1][1:], hist[0]

#
#  Subtact H-alpha background level.
#
#  param[in] hist_bgfg: Histogram of
#  background + foreground level.
#  param[in] hist_bg: Histogram of
#  background level.
#
#  @return Histogram of foreground level.
#
def bg_sub(hist_bgfg, hist_bg):
    
    return hist_bgfg - hist_bg

#
#  Subtact H-alpha background for a
#  single cluster.
#
#  @param[in] gals: List of galaxies.
#  @param[in] clt: DataFrame element of
#  given cluster. 
#  @param[in] opts: List of arguments.
#
#  @return Bins and histogram of
#  foreground level.
#
def single(gals, clt, opts):

    x, hist_bgfg = get_bgfg(gals, clt, opts)
    hist_bg = get_bg(gals, opts)

    return x, (bg_sub(hist_bgfg, hist_bg), )

#
#  Subtact H-alpha background for a
#  stack of clusters.
#
#  @param[in] gals: List of galaxies.
#  @param[in] clt: DataFrame of a sample
#  of clusters. 
#  @param[in] opts: List of arguments.
#  @param[in] m_data: Array of matches.
#
#  @return Bins and histogram of
#  foreground level.
#
def stack(gals, clt, opts, m_data = None):

    count_match = 0
    count_fake = 0
    
    res_match = 0.0
    res_fake = 0.0

    for i in range(clt.id.size):
        
        x, y = get_bgfg(gals, clt.iloc[[i]], opts)

        if np.any(m_data):
            
            if m_data[clt.iloc[[i]].id.base[0]] == 'Y':
                res_match += y
                count_match += 1
            
            else:
                res_fake += y
                count_fake += 1

        else:
            res_match += y
            count_match += 1
    
    hist_bg = get_bg(gals, opts)

    hist_match = None
    hist_fake = None
    
    if count_match >= 1:
        hist_match = bg_sub(res_match / count_match, hist_bg)
    if count_fake >= 1:
        hist_fake = bg_sub(res_fake / count_fake, hist_bg)

    return x, (hist_match, hist_fake), (count_match, count_fake)
