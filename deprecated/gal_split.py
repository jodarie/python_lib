###################
# GAL_SPLIT V.1.0 #
###################

import math
import numpy as np
import numpy.random as npr
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from collections import Counter
import sys, getopt
import os.path


def usage():
   print 'test.py -i <input_file> -n <nsplit> -o <overlap>'
   exit(-1)

def file_name_error(file_name):
      if os.path.isfile(file_name)==False:
         print 'ERROR! Input file name [%s] does not exist.' % file_name
         usage()
      elif file_name=='' or file_name[0][0]=='-':
         print 'ERROR! Input file name not specified.'
         usage()   

def main(argv):

   ##################
   # READ ARGUMENTS #
   ##################

   input_file = ''
   
   try:
      opts, args = getopt.getopt(argv,"i:n:o:h",["ifile=","nsplit=","help"])
   except getopt.GetoptError:
      usage()
   for opt, arg in opts:
      if opt in ('-h','--help'):
         usage()
      elif opt in ("-i", "--ifile"):                                                       #set input file name
         input_file = arg
      elif opt in ("-n", "--nsplit"):                                                      #set number of splits
         nsplit = int(arg)
      elif opt in ("-o", "--overlap"):                                                     #set overlap range
         overlap = float(arg)        
      else:
         print 'Fuck you asshole!'
         exit(0)

   #############
   # READ DATA #
   #############

   print '================================================================='

   file_name_error(input_file)                                                            #make sure files exit
   print 'Reading file:', input_file,        
   file_data = np.genfromtxt(input_file, dtype="S", unpack = True)                        #read file
   print '\tComplete'

   ra = np.array(file_data[1]).astype('float')
   dec = np.array(file_data[2]).astype('float')

   min_ra = min(ra)
   max_ra = max(ra)
   min_dec = min(dec)
   max_dec = max(dec)

   nplit = np.sqrt(float(nsplit))

   print nsplit

   d_ra = max_ra - min_ra
   d_dec = max_dec - min_dec

   print d_ra, d_dec

   print '================================================================='
   
if __name__ == "__main__":
   main(sys.argv[1:])
