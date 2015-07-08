## @file bgsub_main.py
#
#  BGSUB MAIN
#
#  Script for performing background subtraction
#  of H-alpha emitting galaxies.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
import bgsub_bgsub, bgsub_plot

#
#  Calcualte H-alpha foreground level
#  for a single cluster.
#
#  @param[in] clt: List of clusters. 
#  @param[in] gals: List of galaxies.
#  @param[in] opts: List of arguments.
#
#  @return Histogram of background +
#  foreground level.
#
def single(c_data, g_data, opts):
    
    if opts.c_id in c_data.id.values:

        if opts.z_lim[0] < c_data.loc[c_data.id == opts.c_id].z < opts.z_lim[1]:

            hist_data = bgsub_bgsub.single(g_data, c_data.loc[c_data.id == opts.c_id], opts)
            bgsub_plot.plot_nz(hist_data, float(c_data.loc[c_data.id == opts.c_id].z), opts)

        else:

            raise ValueError('Cluster redshift [' + str(float(c_data.loc[c_data.id == opts.c_id].z)) + \
                             '] not within limits ' + str(opts.z_lim) + '.')

    else:
        
        raise ValueError('ID = ' + opts.c_id + ' is not valid.')

#
#  Calcualte H-alpha foreground level
#  for stacked clusters.
#
#  @param[in] clt: List of clusters. 
#  @param[in] gals: List of galaxies.
#  @param[in] opts: List of arguments.
#  @param[in] m_data: List of matches.
#
#  @return Histogram of background +
#  foreground level.
#
def stack(c_data, g_data, opts, m_data = None):

    index = ((c_data.z >= opts.z_lim[0]) &
             (c_data.z < opts.z_lim[1]) &
             (c_data.lamb >= opts.l_lim[0]) &
             (c_data.lamb < opts.l_lim[1]))
    
    hist_data = bgsub_bgsub.stack(g_data, c_data.loc[index], opts, m_data)    
    bgsub_plot.plot_nz(hist_data, 0.0, opts)
