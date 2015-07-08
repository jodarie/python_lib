import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   x, y = np.genfromtxt(inputfile, unpack = True) 
   plt.plot(x,y)
   plt.savefig(outputfile)
   print "Output saved to", outputfile

if __name__ == "__main__":
   main(sys.argv[1:])




