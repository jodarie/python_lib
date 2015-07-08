######################
# CLUSTER_VIEW V.1.0 #
######################

import math
import numpy as np
import numpy.random as npr
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from collections import Counter
import sys, getopt

def usage():
   print 'get it right!'
   sys.exit(2)

def main(argv):
   
   ##################
   # READ ARGUMENTS #
   ##################
   
   try:
      opts, args = getopt.getopt(argv,"i:c:h",["ifile=","cluster=","help"])
   except getopt.GetoptError:
      usage()
   if len(opts) < 1:
      usage()
   for opt, arg in opts:
      if opt in ('-h','--help'):
         usage()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-c", "--cluster"):
         cluster_id = arg
      else:
         usage()

   #############
   # READ DATA #
   #############
         
   data = np.genfromtxt(inputfile, dtype="S", unpack = True)

   cluster_index = np.transpose(np.where(data[0,:] == cluster_id))
   
   gal_ra = np.array(np.transpose(data[4,cluster_index]), dtype='f')[0]   
   gal_dec = np.array(np.transpose(data[5,cluster_index]), dtype='f')[0]
   gal_z = np.array(np.transpose(data[6,cluster_index]), dtype='f')[0]
   halo_id = np.transpose(data[7,cluster_index])[0]

   halo_id_count = Counter(halo_id)

   unique_halo_id = np.unique(halo_id)

   halo_colour = {}
   colour_count = 0
   for val in unique_halo_id:
      halo_colour[val] = colour_count
      colour_count += 1

   halo_colour_val = []
   for val in halo_id:
      halo_colour_val.append(halo_colour[val])

   halo_colour_val = np.array(halo_colour_val)      
   
   ################
   # OUTPUT PLOTS #
   ################

   ra_low = min(gal_ra)
   ra_up = max(gal_ra)
   dec_low = min(gal_dec)
   dec_up = max(gal_dec)
   z_low = min(gal_z)
   z_up = max(gal_z)

   fig1 = plt.figure()
   ax1 = fig1.add_subplot(111, projection='3d')

   fig2 = plt.figure()
   ax2 = fig2.add_subplot(111)

   colours=['b', 'g', 'r', 'c', 'y', 'm', 'indigo', 'olive', 'darkorange', 'brown', 'lightskyblue', 'lime', 'deeppink', 'mediumturquoise', 'gray']
   temp_colour = 0
   for key, val in halo_id_count.iteritems():
      if val > 3 and temp_colour < len(colours):
         index = np.where(halo_id == key)
         temp_ra = np.array(gal_ra[index])
         temp_dec = np.array(gal_dec[index])
         temp_z = np.array(gal_z[index])
         temp_label = "Halo " + key + ' (' + str(val) + ')'
         ax1.plot(temp_ra, temp_dec, temp_z, 'o', c = colours[temp_colour], label = temp_label)
         nz_binsize = 0.005
         nz_nbins = int(math.floor((max(temp_z) - min(temp_z)) / nz_binsize))
         if nz_nbins > 0:
            ax2.hist(temp_z, nz_nbins, normed=1, facecolor = colours[temp_colour], alpha=0.5, label = temp_label)
            temp_colour += 1
      else:
         index = np.where(halo_id == key)
         temp_ra = np.array(gal_ra[index])
         temp_dec = np.array(gal_dec[index])
         temp_z = np.array(gal_z[index])
         ax1.plot(temp_ra, temp_dec, temp_z, 'x', c = 'black')
         
   ax1.set_xlim(ra_low,ra_up)
   ax1.set_ylim(dec_low,dec_up)
   ax1.set_zlim(z_low,z_up)
   ax1.xaxis.set_major_locator(MaxNLocator(5))
   ax1.yaxis.set_major_locator(MaxNLocator(5))
   ax1.zaxis.set_major_locator(MaxNLocator(5))
   ax1.tick_params(axis='both', which='major', labelsize=7)
   ax1.set_title('Cluster Galaxies')
   ax1.set_xlabel('RA')
   ax1.set_ylabel('DEC')
   ax1.set_zlabel('Z')
   ax1.legend(loc='upper left', numpoints=1, ncol=4, fontsize=5, bbox_to_anchor=(0, 0))

   ax2.set_xlabel('$z$')
   ax2.set_title('N(z) Distribution')
   ax2.legend(loc='upper right', numpoints=1, ncol=1, fontsize=5, bbox_to_anchor=(1, 1))

   fig1_name = "cluster_" + cluster_id + "_radecz.pdf"
   fig1.savefig(fig1_name)
   print "Output saved to", fig1_name

   fig2_name = "cluster_" + cluster_id + "_nz.pdf"
   fig2.savefig(fig2_name)
   print "Output saved to", fig2_name

if __name__ == "__main__":
   main(sys.argv[1:])
