##################
# MEM_MATCH v1.0 #
##################

import argparse, os.path
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

import functions
from cluster_class import Clusterx

# READ ARGUMENTS

parser = argparse.ArgumentParser('MEM_MATCH OPTIONS:')
parser.add_argument('-o', '--obs_mem_file', dest = 'obs_mem_file',
                    help = 'Observed members file name.')
parser.add_argument('-m', '--mock_mem_file', dest = 'mock_mem_file',
                    help = 'Mock members file name.')
parser.add_argument('--obs_id', dest = 'obs_id_col', type = int, default = 1, 
                    help = 'Observed cluster ID column number.')
parser.add_argument('--obs_mem_id', dest = 'obs_mem_id_col', type = int, default = 4, 
                    help = 'Observed cluster member ID column number.')
parser.add_argument('--obs_proxy', dest = 'obs_proxy_col', type = int, default = 2, 
                    help = 'Observed cluster mass proxy column number.')
parser.add_argument('--obs_z', dest = 'obs_z_col', type = int, default = 3, 
                    help = 'Observed cluster redshift column number.')
parser.add_argument('--mock_id', dest = 'mock_id_col', type = int, default = 6, 
                    help = 'Mock halo ID column number.')
parser.add_argument('--mock_mem_id', dest = 'mock_mem_id_col', type = int, default = 1, 
                    help = 'Mock halo member ID column number.')
parser.add_argument('--mock_mass', dest = 'mock_mass_col', type = int, default = 7, 
                    help = 'Mock halo mass column number.')
parser.add_argument('--mock_ngal', dest = 'mock_ngal_col', type = int, default = 9, 
                    help = 'Mock halo ngal column number.')
parser.add_argument('--mass_min', dest = 'mass_min', type = float, default = 13.0, 
                    help = 'Minimum mass value for histogram.')
parser.add_argument('--mass_max', dest = 'mass_max', type = float, default = 15.5, 
                    help = 'Maximum mass value for histogram.')
parser.add_argument('--mass_bin', dest = 'mass_bin', type = float, default = 0.1, 
                    help = 'Mass bin size for histogram.')
parser.add_argument('--proxy_min', dest = 'proxy_min', type = float, default = 1.0, 
                    help = 'Minimum mass proxy value for histogram.')
parser.add_argument('--proxy_max', dest = 'proxy_max', type = float, default = 3.0, 
                    help = 'Maximum mass proxy value for histogram.')
parser.add_argument('--proxy_bin', dest = 'proxy_bin', type = float, default = 0.1, 
                    help = 'Mass proxy bin size for histogram.')
parser.add_argument('--z_min', dest = 'z_min', type = float, default = 0.0, 
                    help = 'Minimum redshift value for histogram.')
parser.add_argument('--z_max', dest = 'z_max', type = float, default = 3.0, 
                    help = 'Maximum redshift value for histogram.')
opts = parser.parse_args()

if not opts.obs_mem_file:
    print 'Error: Observed members file not specified.'
    exit()

if not opts.mock_mem_file:
    print 'Error: Mock members file not specified.'
    exit()
 
if not os.path.exists(opts.obs_mem_file):
    print 'Error: Observed members file not found.'
    exit()

if not os.path.exists(opts.mock_mem_file):
    print 'Error: Mock members file not found.'
    exit()

# READ FILES

obs_gals = np.genfromtxt(opts.obs_mem_file, unpack = True, dtype = 'float')
mock_gals = np.genfromtxt(opts.mock_mem_file, unpack = True, dtype = 'float')

obs_mem_id = np.array(obs_gals[opts.obs_mem_id_col - 1], dtype = 'int')
mock_mem_id = np.array(mock_gals[opts.mock_mem_id_col - 1], dtype = 'int')
og_mock_mass = np.array(mock_gals[opts.mock_mass_col - 1], dtype = 'float')

# SORT DATA

index0 = np.argsort(obs_mem_id)

obs_gals = obs_gals[:, index0]
obs_mem_id = obs_mem_id[index0]

index1 = np.argsort(mock_mem_id)

mock_gals = mock_gals[:, index1]
mock_mem_id = mock_mem_id[index1]

# FIND MATCHES

index0 = np.in1d(obs_mem_id, mock_mem_id, assume_unique = True)

obs_gals = obs_gals[:, index0]
obs_mem_id = obs_mem_id[index0]

index1 = np.in1d(mock_mem_id, obs_mem_id, assume_unique = True)

mock_id = mock_gals[opts.mock_id_col - 1, index1]
mock_mass = mock_gals[opts.mock_mass_col - 1, index1]
mock_ngal = mock_gals[opts.mock_ngal_col - 1, index1]

obs_gals = np.vstack([obs_gals, mock_mass, mock_id, mock_ngal])
obs_id = np.unique(np.array(obs_gals[opts.obs_id_col - 1], dtype = 'float'))

# BUILD CLUSTERS

clusters = []
count = 0

for id_val in obs_id:
    clusters.append(Clusterx(count))    
    for member in np.transpose(obs_gals)[(obs_gals[opts.obs_id_col - 1] == id_val)]:
        clusters[count].add_mem(member)
    count += 1

for cluster in clusters:
    cluster.props(opts.obs_proxy_col, opts.obs_z_col)
    cluster.halo_count()
    cluster.mass_hist(opts.mass_min, opts.mass_max, opts.mass_bin)

# HALO MASS HISTOGRAM

hm_hist = functions.get_hist(og_mock_mass, opts.mass_min, opts.mass_max, opts.mass_bin)
    
# MAKE HISTOGRAM MATRIX

n_mass_bins = functions.num_bins(opts.mass_min, opts.mass_max, opts.mass_bin)
n_proxy_bins = functions.num_bins(opts.proxy_min, opts.proxy_max, opts.proxy_bin)

hist_sum = np.zeros(n_mass_bins)
hist_sum_matrix = np.zeros((n_proxy_bins, n_mass_bins))

for cluster in clusters:

    #hist_sum += functions.scale(cluster.hist)
    hist_sum += cluster.hist

    proxy_bin = n_proxy_bins - 1 - \
      functions.find_bin(np.log10(cluster.proxy), opts.proxy_min, opts.proxy_bin)

    for i in range(n_proxy_bins):
        if proxy_bin == i and opts.z_min <= cluster.z < opts.z_max:
            hist_sum_matrix[proxy_bin] += functions.scale(cluster.hist)
                   
for i in range(n_proxy_bins):
    hist_sum_matrix[i] = functions.scale(hist_sum_matrix[i])
        
# SAVE MATRIX TO FILE

output = np.transpose(np.vstack([hm_hist, hist_sum, (hist_sum / hm_hist[1])]))

file_name = opts.obs_mem_file + '.hist.txt'
np.savetxt(file_name, output, fmt = '%.3f')
print 'Data saved to:', file_name

file_name = opts.obs_mem_file + '.matrix.txt'
output = np.fliplr(np.transpose(np.vstack([hist_sum_matrix,
                                           clusters[0].hist_x])))
np.savetxt(file_name, output, fmt = '%.3f')
print 'Data saved to:', file_name

# MAKE PLOTS

z_range = '.' + str(opts.z_min) + '_to_' + str(opts.z_max)

#mp.rcParams.update({'font.size': 6})

plt.figure(1)
extent = [opts.mass_min, opts.mass_max, opts.proxy_min, opts.proxy_max]
plt.imshow(hist_sum_matrix, extent = extent,
           interpolation = 'spline16', cmap = cm.jet,
           aspect = 'auto', norm = LogNorm())
cbar = plt.colorbar()
cbar.set_label(r"$P(\lambda|M)$")
plt.autoscale(False)
plt.xlabel(r"$Log_{10} M$ $[M_{\odot} h^{-1}]$", {'fontsize' : 15})
plt.ylabel(r'$Log_{10} \lambda$', {'fontsize' : 15})
plt_name = opts.obs_mem_file + z_range + '.mass_obs.pdf'
plt.savefig(plt_name)
print "Plot saved to:", plt_name
plt.close()

text = []
for i in range(n_proxy_bins):
    val = opts.proxy_min + (i + 0.5) * opts.proxy_bin
    text.append(r"$\lambda$ = " + str(val))

mp.rcParams.update({'font.size': 10})
        
fig = plt.figure(2)
fig.set_size_inches(30.0, 5.0)
gs = gridspec.GridSpec(1, n_proxy_bins)
for i in range(n_proxy_bins):
    ax = fig.add_subplot(gs[i])
    if i > 0:
        ax.yaxis.set_visible(False)
    else:
        ax.set_ylabel(r"$P(\lambda|M)$")
    ax.set_xlabel(r"$Log_{10} M$ $[M_{\odot} h^{-1}]$")
    ax.set_title(text[i])
    ax.plot(clusters[0].hist_x, np.flipud(hist_sum_matrix)[i], '-')
    ax.set_xlim(opts.mass_min, opts.mass_max)
    ax.set_ylim(0.0, 1.0)
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(3))
plt_name = opts.obs_mem_file + z_range + '.p_lambda.pdf'
plt.savefig(plt_name)
print "Plot saved to:", plt_name
plt.close()
