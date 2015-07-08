#! /Users/Rowen/Documents/Library/anaconda/bin/python

################
# SPEARMAN RHO #
################

import math
import numpy as np
import scipy.stats as ss
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Read Data

AG = np.genfromtxt('AG_peak', dtype="float", unpack = True)
CB = np.genfromtxt('CB', dtype="float", unpack = True)
FB = np.genfromtxt('FB2', dtype="float", unpack = True)
FD = np.genfromtxt('FD', dtype="float", unpack = True)
RL = np.genfromtxt('RL', dtype="float", unpack = True)
SF = np.genfromtxt('SF', dtype="float", unpack = True)

min_z = 0.1
max_z = 2.6
bin_z = 0.5

n_z_bins = int(math.floor((max_z - min_z) / bin_z))
x_z_vals = (np.arange(n_z_bins) + 0.5) * bin_z + min_z

# Spearman

def find_spear(data, n_z_bins, min_z, bin_z):
    z = data[0]
    mass = data[1]
    proxy = data[2]    
    z_bin_index = np.floor(np.round((z - min_z) / bin_z, 8)).astype('int')
    for i in range(n_z_bins):
        m_index = mass[z_bin_index == i]
        p_index = proxy[z_bin_index == i]            
        if m_index.size > 1:
            rho, p1 = ss.spearmanr(m_index, p_index)
            rho_err = 0.6325 / (len(m_index) - 1) ** 0.5
        else:
            rho = 0.0
            rho_err = 0.0
        if i == 0:
            rhos = np.array(rho)
            rhos_err = np.array(rho_err)
        else:
            rhos = np.hstack((rhos, rho))
            rhos_err = np.hstack((rhos_err, rho_err))
    return rhos, rhos_err

AG_rho, AG_rho_err = find_spear(AG, n_z_bins, min_z, bin_z)
CB_rho, CB_rho_err = find_spear(CB, n_z_bins, min_z, bin_z)
FB_rho, FB_rho_err = find_spear(FB, n_z_bins, min_z, bin_z)
FD_rho, FD_rho_err = find_spear(FD, n_z_bins, min_z, bin_z)
RL_rho, RL_rho_err = find_spear(RL, n_z_bins, min_z, bin_z)
SF_rho, SF_rho_err = find_spear(SF, n_z_bins, min_z, bin_z)

#z vs. spearman rho
 
fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.plot(x_z_vals, AG_rho, linewidth=1.0, c='green', label='AG')
ax.plot(x_z_vals, CB_rho, linewidth=1.0, c='#FFBF00', label='CB')
ax.plot(x_z_vals, FB_rho, linewidth=1.0, c='blue', label='FB')
ax.plot(x_z_vals, FD_rho, linewidth=1.0, c='magenta', label='FD')
ax.plot(x_z_vals, RL_rho, linewidth=1.0, c='purple', label='RL')
ax.plot(x_z_vals, SF_rho, linewidth=1.0, c='red', label='SF')
ax.errorbar(x_z_vals, AG_rho, yerr=AG_rho_err, fmt='green')
ax.errorbar(x_z_vals, CB_rho, yerr=CB_rho_err, fmt='#FFBF00')
ax.errorbar(x_z_vals, FB_rho, yerr=FB_rho_err, fmt='blue')
ax.errorbar(x_z_vals, FD_rho, yerr=FD_rho_err, fmt='magenta')
ax.errorbar(x_z_vals, RL_rho, yerr=RL_rho_err, fmt='purple')
ax.errorbar(x_z_vals, SF_rho, yerr=SF_rho_err, fmt='red')
ax.set_xlim(min_z, max_z)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel('z')
ax.set_ylabel(r'$\rho$')
ax.set_title('Spearman Rank Order Coeffcient')
ax.legend()
fig_name = 'test_spearman.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name
