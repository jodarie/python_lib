#########################
# IDL2ASCII V.1.0       #
# Samuel Farrens (2013) #
#########################

import os.path
import sys, getopt
import numpy as np
from scipy.io.idl import readsav

def usage():
   print 'FOR HELP:'
   print '>> python -h'
   print 'BASIC USAGE:'
   print '>> python idl2ascii -i <input_file> -o <output_file>'
   print '(Note: This script requires numpy and scipy)'
   sys.exit(2)

def help():
   print '|-------------------------------------------------------------------------------------|'
   print '| This code converts a given IDL SAVE (.sav) file to an ASCII table.                  |'
   print '|-------------------------------------------------------------------------------------|'  
   print '| OPTIONS:                                                                            |'
   print '|                                                                                     |'
   print '| -h [--help]    Displays this help page.                                             |'
   print '|                                                                                     |'
   print '| -i [--input]   An input file name should be provided as an argument following this  |'
   print '|                option.                                                              |'
   print '|                                                                                     |'
   print '| -o [--output]  An output file name can be provided as an argument following this    |'
   print '|                option.                                                              |'
   print '|                                                                                     |'
   print '| -c [--col]     The IDL column names can be specified following this option.         |'
   print '|                                                                                     |'
   print '| -v [--verbose] This option displays the IDL SAVE file content verbose.              |'
   print '|                                                                                     |'
   print '|-------------------------------------------------------------------------------------|'
   print '| EXAMPLES:                                                                           |'
   print '|                                                                                     |'
   print '|  > python idl2ascii -i <input_file_name> -v                                         |'
   print '|                                                                                     |'
   print '| This will read an an IDL SAVE file and output the content verbose.                  |'
   print '|                                                                                     |'
   print '|  > python idl2ascii -i <input_file_name> -o <output_file_name>                      |'
   print '|                                                                                     |'
   print '| This will read an an IDL SAVE file and output all of the content to an ASCII file.  |'
   print '|                                                                                     |'
   print '|  > python idl2ascii -i <input_file_name> -o <output_file_name> -c <column_name1>    |'
   print '|  ... -c <column_name2> -c <column_name3>                                            |'
   print '|                                                                                     |'
   print '| This will read an an IDL SAVE file and output only the specified content to an      |'
   print '| ASCII file.                                                                         |'
   print '|                                                                                     |'
   print '|-------------------------------------------------------------------------------------|'
   sys.exit(2)

def main(argv):
   
   #####################
   # READ ARGUMENTS    #
   #####################

   inputfile=''
   outputfile=''
   columns=[]
   
   try:
      opts, args = getopt.getopt(argv,"i:o:c:hv",["input=","output=","help","verbose"])
   except getopt.GetoptError:
      usage()
   if len(opts)<1:
      usage()
   verb = False
   for opt, arg in opts:
      if opt in ('-h','--help'):
         help()
      elif opt in ('-i', '--input'):
         inputfile = arg
      elif opt in ('-o', '--output'):
         outputfile = arg
      elif opt in ('-c', '--col'):
         columns.append(arg)
      elif opt in ('-v','--verbose'):
         verb = True
         
   #####################
   # READ IDL FILE     #
   #####################

   if os.path.isfile(inputfile)==False:
      print 'ERROR! Input file name [%s] does not exist.' % inputfile
      usage()
         
   if inputfile=='' or inputfile[0][0]=='-':
      print 'ERROR! Input file name not specified.'
      usage()
      
   data = readsav(inputfile,verbose=verb)
   key = data.keys()

   #####################
   # WRITE ASCII FILE  #
   #####################

   if outputfile=='':
      sys.exit(1)
   
   out = open(outputfile,'w')

   match=0
   if len(columns)>0:
      columns = np.array(columns)
      cols = np.empty(len(columns),dtype=int)
      print>> out,'#',
      for i in range(len(columns)):
         for j in range(len(key)):
            if(key[j]==columns[i]):
               print>>out,key[j],
               cols[i]=j
               match+=1
      print>> out,''
      if(match!=len(cols)):
         print 'ERROR! One or more column names not found in file.'
         usage()
      for i in range(len(data[key[0]])):
         for j in range(len(cols)):
            print>> out,data[key[cols[j]]][i],
         print>> out,''

   else:
      print>> out,'#',
      for i in range(len(key)):
         print>>out,key[i],
      print>> out,''
      for i in range(len(data[key[0]])):
         for j in range(len(key)):
            print>> out,data[key[j]][i],
         print>> out,''
      
if __name__ == "__main__":
    main(sys.argv[1:])
