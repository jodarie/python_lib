## @file bgsub_plot.py
#
#  BGSUB PLOT FUNCTIONS 
#
#  Functions for plotting.
#  
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
import matplotlib.pyplot as plt

##
#  Plot N(z) for background 
#  subtracted H-alpha galaxies.
#
#  @param[in] hist_data: Histogram data.
#  @param[in] z: Cluster redshift.
#  @param[in] opts: List of arguments.
#
#  @return DataFrame of clusters.
#
def plot_nz(hist_data, z, opts):
        
    if len(hist_data) > 2:
        label1 = r'FG H$_{\alpha}$ Matched: ' + str(hist_data[2][0])
        label2 = r'FG H$_{\alpha}$ Unmatched: ' + str(hist_data[2][1])
        label3 = r'z$_c$ - z = %.3f' % z
    else:
        label1 = r'FG H$_{\alpha}$'
        label3 = 'z = %.3f' % z
        
    plt.step(hist_data[0], hist_data[1][0], 'r', linestyle = '-', label = label1)

    if len(hist_data[1]) > 1 and np.any(hist_data[1][1]):
        plt.step(hist_data[0], hist_data[1][1], 'b', linestyle = '--', label = label2)
        limits = (min(np.nanmin(hist_data[1][0]), np.nanmin(hist_data[1][1])) - 1,
              max(np.nanmax(hist_data[1][0]), np.nanmax(hist_data[1][1])) + 1)
    else:
        limits = (np.nanmin(hist_data[1][0]) - 1, np.nanmax(hist_data[1][0]) + 1)

    plt.plot((z, z), limits, 'k:', label = label3)
    
    plt.legend(prop = {'size':10})
        
    plt.ylim(*limits)

    if opts.c_id:
        file_name = opts.input_file[1] + '_cluster_%s_nz.pdf' % opts.c_id
        plt.xlabel('z', {'fontsize' : 20})
        plt.ylabel('N(z)', {'fontsize' : 20})
        plt.title('Cluster Candidate ID:%s' % opts.c_id)
        plt.xlim(*opts.z_lim)
    else:
        file_name = opts.input_file[1] + '_clusters_lambda.%.2f_nz.pdf' % opts.l_lim[0]
        plt.xlabel('z - z$_c$', {'fontsize' : 20})
        plt.ylabel('N(z - z$_c$)', {'fontsize' : 20})
        plt.title(r'Range: ' + str(opts.l_lim[0]) + '$\leq\lambda\leq$' + str(opts.l_lim[1]))
        
    plt.savefig(file_name)
    print 'Plot saved to:', file_name
