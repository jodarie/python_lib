################
# PYDIAG V.2.0 #
################

import sys, getopt, math, os.path
import numpy as np, numpy.random as npr, matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from collections import Counter
import matplotlib.patches as patches


def usage():
   '''
   Function that prints basic code usage.
   '''
   print 'To run:'
   print '>> pydiag -p pydiag_param.ini'
   print 'Where pydiag_param.ini is a file containing the pydiag parameters parameter.'
   exit(-1)

def file_name_error(file_name):
   '''
   Function that returns an error message for missing files
   '''
   if os.path.isfile(file_name)==False:
      print 'ERROR! Input file name [%s] does not exist.' % file_name
      usage()
   elif file_name=='' or file_name[0][0]=='-':
      print 'ERROR! Input file name not specified.'
      usage()   

def get_rich(galaxy_list,unique_galaxy_list):
   '''
   Function to determine richness of a given halo/cluster
   '''
   galaxy_counter = Counter(galaxy_list)
   rich_list = []
   for i in range(len(unique_galaxy_list)):
      rich_list.append(galaxy_counter[unique_galaxy_list[i]])
   rich_list = np.array(rich_list)
   return rich_list
   
def find_galaxy_match(cluster_galaxy_ids,mock_galaxy_ids):
   '''
   Function to match galaxy IDs between cluster files and mocks
   '''
   g_tf = np.in1d(mock_galaxy_ids,cluster_galaxy_ids)      
   g_index = np.arange(mock_galaxy_ids.shape[0])[g_tf]
   return g_index

def find_halo_match(halo_ids,cluster_ids,cluster_now):
   '''
   Function to find the most common Halo ID associated to a given cluster
   '''
   c_index = cluster_ids==cluster_now
   c_temp = halo_ids[c_index]
   c_count = Counter(c_temp)
   c_common = c_count.most_common(1)
   return c_common[0][0],c_common[0][1]

def loop_through_matches(halo_ids,cluster_ids,unique_cluster_ids):
   '''
   Function that loops through all cluster to halo matches
   '''
   m_id_match = []
   m_num_match = []
   for i in range(len(unique_cluster_ids)):                                               
      m_id_match_now,m_num_match_now = find_halo_match(halo_ids,cluster_ids,unique_cluster_ids[i])
      m_id_match.append(m_id_match_now)
      m_num_match.append(m_num_match_now)
   m_id_match = np.array(m_id_match,dtype='S')
   m_num_match = np.array(m_num_match,dtype='i')
   return m_id_match,m_num_match

def nan2one(array):
   '''
   Convert nan values to 1.0
   '''
   array[np.isnan(array)] = 1.0

def main(argv):

   print '================================================================='

   ###################
   # READ PARAMETERS #
   ###################

   def read_params(param_file):
      '''
      Function to read parameters from the parameter file.
      '''
      print 'Reading file:', param_file,        
      parameter_data = np.genfromtxt(param_file, dtype="S", unpack = True)                 
      print '\t\t\tComplete'
      if 'cluster_file' in parameter_data:
         index = int(np.where(parameter_data == 'cluster_file')[1])
         cluster_file = parameter_data[1,index]
      else: cluster_file = ''
      if 'mock_file' in parameter_data:
         index = int(np.where(parameter_data == 'mock_file')[1])
         mock_file = parameter_data[1,index]
      else: mock_file = ''
      if 'c_id_col' in parameter_data:
         index = int(np.where(parameter_data == 'c_id_col')[1])
         c_id_col = int(parameter_data[1,index]) - 1
      else: c_id_col = -2
      if 'cg_id_col' in parameter_data:
         index = int(np.where(parameter_data == 'cg_id_col')[1])
         cg_id_col = int(parameter_data[1,index]) - 1
      if 'cg_z_col' in parameter_data:
         index = int(np.where(parameter_data == 'cg_z_col')[1])
         cg_z_col = int(parameter_data[1,index]) - 1   
      else: cg_id_col = -2
      if 'm_id_col' in parameter_data:
         index = int(np.where(parameter_data == 'm_id_col')[1])
         m_id_col = int(parameter_data[1,index]) - 1
      else: m_id_col = -2
      if 'mg_id_col' in parameter_data:
         index = int(np.where(parameter_data == 'mg_id_col')[1])
         mg_id_col = int(parameter_data[1,index]) - 1
      if 'mg_z_col' in parameter_data:
         index = int(np.where(parameter_data == 'mg_z_col')[1])
         mg_z_col = int(parameter_data[1,index]) - 1   
      else: mg_id_col = -2
      if 'm_mass_col' in parameter_data:
         index = int(np.where(parameter_data == 'm_mass_col')[1])
         m_mass_col = int(parameter_data[1,index]) - 1
      else: m_mass_col = -2
      if 'match_threshold' in parameter_data:
         index = int(np.where(parameter_data == 'match_threshold')[1])
         match_threshold = float(parameter_data[1,index])
      else: match_threshold = -2
      if 'rich_bin_size' in parameter_data:
         index = int(np.where(parameter_data == 'rich_bin_size')[1])
         rich_bin_size = float(parameter_data[1,index])
      else: rich_bin_size = -2
      if 'max_rich_bin' in parameter_data:
         index = int(np.where(parameter_data == 'max_rich_bin')[1])
         max_rich_bin = float(parameter_data[1,index])
      else: max_rich_bin = -2
      if 'min_rich_bin' in parameter_data:
         index = int(np.where(parameter_data == 'min_rich_bin')[1])
         min_rich_bin = float(parameter_data[1,index]) 
      else: min_rich_bin = -2
      if 'z_bin_size' in parameter_data:
         index = int(np.where(parameter_data == 'z_bin_size')[1])
         z_bin_size = float(parameter_data[1,index]) 
      else: z_bin_size = -2
      if 'max_z_bin' in parameter_data:
         index = int(np.where(parameter_data == 'max_z_bin')[1])
         max_z_bin = float(parameter_data[1,index])
      else: max_z_bin = -2
      if 'min_z_bin' in parameter_data:
         index = int(np.where(parameter_data == 'min_z_bin')[1])
         min_z_bin = float(parameter_data[1,index])
      else: min_z_bin = -2
      return cluster_file,mock_file,c_id_col,cg_id_col,cg_z_col,m_id_col,mg_id_col,mg_z_col,m_mass_col,match_threshold,rich_bin_size,max_rich_bin,min_rich_bin,z_bin_size,max_z_bin,min_z_bin
   
   ##################
   # READ ARGUMENTS #
   ##################

   cluster_file = ''
   mock_file = ''
   c_id_col = -1
   cg_id_col = -1
   cg_z_col = -1
   m_id_col = -1
   mg_id_col = -1
   mg_z_col = -1
   m_mass_col = -1
   match_threshold = -1
   rich_bin_size = -1
   max_rich_bin = -1
   min_rich_bin = -1
   z_bin_size = -1
   max_z_bin = -1
   min_z_bin = -1
   
   try:
      opts, args = getopt.getopt(argv,"c:m:p:a:b:c:d:e:f:g:h:i:j:k:l:z:o:p",["cfile=","hfile=","pfile=","c_id_col=","cg_id_col=","cg_z_col="
                                                               "m_id_col=","mg_id_col=","mg_z_col=","m_mass_col=","binsize=",
         "match_threshold=","rich_bin=","max_rich_bin=","min_rich_bin=","z_bin_size=","max_z_bin","min_z_bin","help"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         usage()
         sys.exit(2)
      elif opt in ("-p", "--pfile"):                                                       #set parameter file name
         param_file = arg
         cluster_file,mock_file,c_id_col,cg_id_col,cg_z_col,m_id_col,mg_id_col,mg_z_col,m_mass_col,match_threshold,rich_bin_size,max_rich_bin,min_rich_bin,z_bin_size,max_z_bin,min_z_bin = read_params(param_file)
      elif opt in ("-c", "--cfile"):                                                       #set cluster galaxy file name
         cluster_file = arg
      elif opt in ("-m", "--mfile"):                                                       #set mock galaxy file name
         mock_file = arg 
      elif opt in ("-a", "--c_id_col"):                                                    #set cluster ID column number
         c_id_col = int(arg) - 1
      elif opt in ("-b", "--cg_id_col"):                                                   #set cluster galaxy ID column number
         cg_id_col = int(arg) - 1
      elif opt in ("-k", "--cg_z_col"):                                                    #set cluster galaxy redshift column number
         cg_z_col = int(arg) - 1
      elif opt in ("-c", "--m_id_col"):                                                    #set halo ID column number
         m_id_col = int(arg) - 1
      elif opt in ("-d", "--mg_id_col"):                                                   #set mock galaxy ID column number
         mg_id_col = int(arg) - 1
      elif opt in ("-l", "--mg_z_col"):                                                    #set mock galaxy redshift column number
         mg_z_col = int(arg) - 1 
      elif opt in ("-e", "--m_mass_col"):                                                  #set halo mass column number
         m_mass_col = int(arg) - 1
      elif opt in ("-f", "--match_threshold"):                                             #set threshold for cluster to halo match
         match_threshold = float(arg)
      elif opt in ("-g", "--rich_bin"):                                                    #set richness bin size
         rich_bin_size = float(arg)
      elif opt in ("-i", "--max_rich_bin"):                                                #set maximum richness
         max_rich_bin = float(arg)
      elif opt in ("-j", "--min_rich_bin"):                                                #set minimum richness
         min_rich_bin = float(arg)
      elif opt in ("-z", "--z_bin_size"):                                                  #set redshift bin size
         z_bin_size = float(arg)
      elif opt in ("-o", "--max_z_bin"):                                                   #set maximum redshift
         max_z_bin = float(arg)
      elif opt in ("-p", "--min_z_bin"):                                                   #set minimum redshift
         min_z_bin = float(arg)
      else:
         print 'Something is wrong!'
         sys.exit(2)

   if match_threshold < 0: match_threshold = 50
   elif match_threshold <= 1: match_threshold = match_threshold * 100

   #############
   # READ DATA #
   #############

   file_name_error(cluster_file)                                                          #make sure files exit
   file_name_error(mock_file)
   print 'Reading file:', cluster_file,        
   cluster_data = np.genfromtxt(cluster_file, dtype="S", unpack = True)                   #cluster galaxies file
   print '\tComplete'
   print 'Reading file:', mock_file,        
   mock_data = np.genfromtxt(mock_file, dtype="S", unpack = True)                         #mock halo galaxies file
   print '\tComplete'

   ################
   # FIND MATCHES #
   ################   

   #from cluster galaxies file
   cg_id = np.array(cluster_data[cg_id_col, :], dtype='i')                                #galaxy id in cluster file   
   cg_index = np.argsort(cg_id)                                                           #sort galaxy ids
   cg_id = cg_id[cg_index]
   c_id = cluster_data[c_id_col, cg_index]                                                #cluster ids
   cg_z = np.array(cluster_data[cg_z_col, cg_index], dtype='f')
   c_id_u, c_id_u_index = np.unique(c_id, return_index=True)                              #unique list of cluster ids   
   c_rich_u = get_rich(c_id, c_id_u)                                                      #cluster richnesses
   c_z_u = cg_z[c_id_u_index]                                                             #cluster redshift
         
   #from mock halo galaxies file
   mg_id = np.array(mock_data[mg_id_col,:],dtype='i')                                     #galaxy id in mock file   
   mg_index = np.argsort(mg_id)                                                           #sort galaxy ids
   mg_id = mg_id[mg_index]
   m_id = mock_data[m_id_col,mg_index]                                                    #halo ids
   if m_mass_col>-1: m_mass = mock_data[m_mass_col, mg_index]                             #halo masses
   m_id_counter = Counter(m_id)                                                           #halo id counter (gives halo richnesses)
   m_id_u, m_id_u_index = np.unique(m_id,return_index=True)                               #unique list of halo ids
   m_rich_u = get_rich(m_id, m_id_u)
   if min_rich_bin < 0: min_rich_bin = min(c_rich_u)
   m_id_u_3 = m_id_u[(m_rich_u >= min_rich_bin)]
   m_rich_u_3 = get_rich(m_id, m_id_u_3)
   
   mg_z = np.array(mock_data[mg_z_col, m_id_u_index], dtype='f')
   m_z_u_3 = mg_z[(m_rich_u >= min_rich_bin)]

   m_in_c = find_galaxy_match(cg_id,mg_id)                                                #index of halos matched to galaxies in cluster file
   m_id_c = m_id[m_in_c]                                                                  #halo ids matched to galaxies in cluster file

   if len(m_id_c) < 1:
      print 'ERROR! No galaxy IDs matched!'
      exit(-2)
      
   m_id_match,m_num_match = loop_through_matches(m_id_c,c_id,c_id_u)                      #find halo match for each galaxy cluster

   num_non_unique = Counter(m_id_match)                                                   #count number of matched haloes
   most_common = np.transpose(num_non_unique.most_common())                               #find number of non-uniquely matched haloes
   non_unique_index = most_common[1] != '1'
   m_id_non_unique = most_common[0][non_unique_index]                                     #halo ids of non-uniquely matched haloes
   m_id_non_unique_flag = np.zeros(len(m_id_non_unique))                                  #flag highest number of matching galaxies
   for i in range(len(m_id_match)):
      for j in range(len(m_id_non_unique)):
         if m_id_match[i] == m_id_non_unique[j]:
            if m_id_non_unique_flag[j] < m_num_match[i]:
               m_id_non_unique_flag[j] = m_num_match[i]
         
   c_one_way = np.zeros(len(c_id_u))                                                      #flag cluster if it satisfies 1-way macth
   c_two_way = np.zeros(len(c_id_u))                                                      #flag cluster if it satisfies 2-way macth
   m_one_way = np.zeros(len(m_id_u_3))                                                    #flag halo if it satisfies 1-way macth
   m_two_way = np.zeros(len(m_id_u_3))                                                    #flag halo if it satisfies 2-way macth
      
   for i in range(len(c_id_u)):                                                           #count number of matches
      m_rich_now = int(m_id_counter[m_id_match[i]])                                       #find matched halo richness
      m_complete = (float(m_num_match[i])/m_rich_now)*100                                 #find the retreived completeness of the matched halo
      c_pure = (float(m_num_match[i])/float(c_rich_u[i]))*100                             #find the purity of the current cluster
      if m_rich_now > 2 and m_complete >= match_threshold:                                #make sure halo has at least 3 members and is complete
         if m_id_match[i] in m_id_non_unique:                                             #check if the cluster is uniquely matched to a halo
            if m_num_match[i] in m_id_non_unique_flag:                                    #if not check if this is the best match for the halo
               c_one_way[i] = 1
               if c_pure >= match_threshold: c_two_way[i] = 1                             #check if the current cluster is pure enough
         else:
            c_one_way[i] = 1
            if c_pure >= match_threshold: c_two_way[i] = 1
               
   if rich_bin_size < 1: rich_bin_size = 1
   if z_bin_size < 1: z_bin_size = 0.05
   if max_z_bin < 1: max_z_bin = 1.2
   if min_z_bin < 1: min_z_bin = 0.0
      
      #if max_rich_bin < 0: max_rich_bin = max(m_rich_u_3)
      #if max_rich_bin < 0: max_rich_bin = max(max(c_rich_u),max(m_rich_u_3))
      
   n_c_rich_bins = int(math.floor((max(c_rich_u) - min(c_rich_u)) / rich_bin_size)) + 1       #set number of richness bins
   n_m_rich_bins = int(math.floor((max(m_rich_u_3) - min(m_rich_u_3)) / rich_bin_size)) + 1   #set number of richness bins
   n_z_bins = int(math.floor((max_z_bin - min_z_bin) / z_bin_size)) + int(math.ceil(((max_z_bin - min_z_bin) / z_bin_size) % 1))
                                                                                              #set number of z bins

   x_c_rich_vals = min(c_rich_u) + rich_bin_size * (np.linspace(0, n_c_rich_bins - 1, n_c_rich_bins) + 0.5)
   x_m_rich_vals = min(m_rich_u_3) + rich_bin_size * (np.linspace(0, n_m_rich_bins - 1, n_m_rich_bins) + 0.5)
   y_z_vals = np.arange(min_z_bin, max_z_bin + z_bin_size, z_bin_size)

   c_bin_count = np.zeros((n_c_rich_bins, n_z_bins)).astype('int')
   m_bin_count = np.zeros((n_m_rich_bins, n_z_bins)).astype('int')
   c_match_bin_count_1 = np.zeros((n_c_rich_bins, n_z_bins)).astype('int')
   c_match_bin_count_2 = np.zeros((n_c_rich_bins, n_z_bins)).astype('int')
   m_match_bin_count_1 = np.zeros((n_m_rich_bins, n_z_bins)).astype('int')
   m_match_bin_count_2 = np.zeros((n_m_rich_bins, n_z_bins)).astype('int')
   
   c_rich_bin_index = np.floor((c_rich_u - min_rich_bin) / rich_bin_size).astype('int')
   m_rich_bin_index = np.floor((m_rich_u_3 - min_rich_bin) / rich_bin_size).astype('int')
   
   c_z_bin_index = np.floor((c_z_u - min_z_bin) / z_bin_size).astype('int')
   m_z_bin_index = np.floor((m_z_u_3 - min_z_bin) / z_bin_size).astype('int')

   for i in range(len(c_id_u)):                                                          #cluster match count in richness bins
      c_bin_count[c_rich_bin_index[i], c_z_bin_index[i]] += 1
      if c_one_way[i] == 1:
         c_match_bin_count_1[c_rich_bin_index[i], c_z_bin_index[i]] += 1
      if c_two_way[i] == 1:
         c_match_bin_count_2[c_rich_bin_index[i], c_z_bin_index[i]] += 1

   for i in range(len(m_id_u_3)):                                                        #halo match count in richness bins
      m_bin_count[m_rich_bin_index[i], m_z_bin_index[i]] += 1
      if m_id_u_3[i] in m_id_match:
         an_index = np.transpose(np.where(m_id_match == m_id_u_3[i]))
         for j in range(len(an_index)):
            if c_one_way[an_index[j]] == 1:
               m_match_bin_count_1[m_rich_bin_index[i], m_z_bin_index[i]] += 1
            if c_two_way[an_index[j]] == 1:
               m_match_bin_count_2[m_rich_bin_index[i], m_z_bin_index[i]] += 1    

   pure_bin_1 = np.array(c_match_bin_count_1).astype("float") / np.array(c_bin_count).astype("float")
   nan2one(pure_bin_1)
   complete_bin_1 = np.array(m_match_bin_count_1).astype("float") / np.array(m_bin_count).astype("float")
   nan2one(complete_bin_1)
   pure_bin_2 = np.array(c_match_bin_count_2).astype("float") / np.array(c_bin_count).astype("float")
   nan2one(pure_bin_2)
   complete_bin_2 = np.array(m_match_bin_count_2).astype("float") / np.array(m_bin_count).astype("float")
   nan2one(complete_bin_2)

   #g1_bin = np.array(1.0 - np.sqrt(((1-pure_bin_2)**2.0 + (1.0 - complete_bin_2)**2.0) / 2.0)).astype("float")
   #nan2one(g1_bin)
   #g2_bin = np.array((complete_bin_2 / complete_bin_1) * (pure_bin_2 / pure_bin_1)).astype("float")
   #nan2one(g2_bin)

   ####################STOP!!!!#######################  
                             
   ##################
   # OUTPUT MATCHES #
   ##################

   match_out_file = cluster_file + '_matching.txt'                                           #file with list of cluster to halo matches
   match_out = open(match_out_file,'w')
   
   if m_mass_col>-1:
      print>> match_out, '#c_id \t\tc_rich \thalo_id \t\thalo_rich \thalo_mass \tmatch_count \thalo_% \tcluster_% \t1-way_match \t2-way_match'
   else:
      print>> match_out, '#c_id \t\tc_rich \thalo_id \t\thalo_rich \tmatch_count \thalo_% \tcluster_% \t1-way_match \t2-way_match'
   for i in range(len(c_id_u)):
      m_rich_now = int(m_id_counter[m_id_match[i]])
      if m_mass_col>-1: m_mass_now = float(m_mass[np.where(m_id==m_id_match[i])[0][0]])
      m_complete = (float(m_num_match[i])/m_rich_now)*100
      c_pure = (float(m_num_match[i])/float(c_rich_u[i]))*100      
      print>> match_out, '%(c_id_u_val)12s \t%(c_rich_u_val)06d \t' % {'c_id_u_val': c_id_u[i], 'c_rich_u_val': c_rich_u[i]},
      print>> match_out, '%(m_id_match_val)s \t%(m_rich_now_val)06d \t\t' % {'m_id_match_val': m_id_match[i], 'm_rich_now_val': m_rich_now},
      if m_mass_col>-1: print>> match_out, '%(m_mass_now_val)6.3f \t\t' % {'m_mass_now_val': m_mass_now},
      print>> match_out, '%(m_num_match_val)06d \t\t' % {'m_num_match_val': m_num_match[i]},
      print>> match_out, '%(m_complete_val)5.1f \t%(c_pure_val)5.1f \t\t' % {'m_complete_val': m_complete, 'c_pure_val': c_pure},
      print>> match_out, '%(one_way_val)3.1f \t\t%(two_way_val)3.1f \t' % {'one_way_val': c_one_way[i], 'two_way_val': c_two_way[i]}

      #match_bin_file = cluster_file + '_matching_bin.txt'                                      #cluster metrics as a function of richness
      # match_bin_out = open(match_bin_file,'w')

      #print>> match_bin_out, '#rich \tc_count  m_count \tcmbc1 \t cmbc2 \t\tmmbc1 \t mmbc2 \t\tcb1 \t cb2 \tpb1 \t pb2 \tg1 \t g2'
      #for i in range(n_rich_bins):
      #print>> match_bin_out, '%(x_rich)04d \t' % {'x_rich': x_rich_vals[i]},
      # print>> match_bin_out, '%(cbc)06d \t %(mbc)06d \t' % {'cbc': c_bin_count[i], 'mbc': m_bin_count[i]},
      # print>> match_bin_out, '%(cmbc1)06d \t %(cmbc2)06d \t' % {'cmbc1': c_match_bin_count_1[i], 'cmbc2': c_match_bin_count_2[i]},
      #print>> match_bin_out, '%(mmbc1)06d \t %(mmbc2)06d \t' % {'mmbc1': m_match_bin_count_1[i], 'mmbc2': m_match_bin_count_2[i]}, 
      #print>> match_bin_out, '%(cb1)5.3f \t %(cb2)5.3f \t' % {'cb1': complete_bin_1[i], 'cb2': complete_bin_2[i]},
      #print>> match_bin_out, '%(pb1)5.3f \t %(pb2)5.3f \t' % {'pb1': pure_bin_1[i], 'pb2': pure_bin_2[i]},
      #print>> match_bin_out, '%(g1)5.3f \t %(g2)5.3f' % {'g1': g1_bin[i], 'g2': g2_bin[i]}
      
   #################
   # PRINT RESULTS #
   #################

   completeness_1 = c_one_way.sum() / len(m_id_u_3)                                         #one-way completeness
   purity_1 = c_one_way.sum() / len(c_id_u)                                                 #one-way purity
   completeness_2 = c_two_way.sum() / len(m_id_u_3)                                         #two-way completeness
   purity_2 = c_two_way.sum() / len(c_id_u)                                                 #two-way purity
   g1_parameter = 1 - np.sqrt(((1-purity_2)**2 + (1-completeness_2)**2) / 2)                #balance between completeness and purity
   g2_parameter = (completeness_2 / completeness_1) * (purity_2 / purity_1)                 #level of merging/fragmentation

   print '================================================================='
   print 'Total Number of Clusters:\t\t\t%07d' % len(c_id_u)
   print 'Total Number of Mock Haloes:\t\t\t%07d' % len(m_id_u)
   print 'Total Number of Mock Haloes (Ngal >= 3):\t%07d' % len(m_id_u_3)
   print '-----------------------------------------------------------------'
   print 'Total Number of One-Way Matches:\t\t%07d' % int(c_one_way.sum())
   print 'One-Way Completeness:\t\t\t\t%5.5f' % completeness_1
   print 'One-Way Purity:\t\t\t\t\t%5.5f' % purity_1
   print '-----------------------------------------------------------------'
   print 'Total Number of Two-Way Matches:\t\t%07d' % int(c_two_way.sum())
   print 'Two-Way Completeness:\t\t\t\t%5.5f' % completeness_2
   print 'Two-Way Purity:\t\t\t\t\t%5.5f' % purity_2
   print '-----------------------------------------------------------------'
   print 'G1 Parameter:\t\t\t\t\t%5.5f' % g1_parameter
   print 'G2 Parameter:\t\t\t\t\t%5.5f' % g2_parameter
   print '================================================================='

   ################
   # OUTPUT PLOTS #
   ################

   xyline = (0,1000)

   for i in range(n_z_bins):

      fig1 = plt.figure() 
      plt.rc('legend',**{'fontsize':6})
      gs1 = gridspec.GridSpec(1, 1)
      ax1 = fig1.add_subplot(gs1[0])

      majorFormatter = FormatStrFormatter('%4.2f')
      ax1.yaxis.set_major_formatter(majorFormatter)
      ax1.plot(xyline, (1, 1), '--', c='0.5', linewidth=1.0)
      ax1.plot(x_m_rich_vals, complete_bin_1[:, i], linewidth=1.5, c='b', label='One-Way')                #one-way completeness plot
      ax1.plot(x_m_rich_vals, complete_bin_2[:, i], linewidth=1.5, c='r', label='Two-Way')                #two-way completeness plot
      if (np.where(m_bin_count[:, i] > 0)[0]).size > 0:
         ax1.plot(x_m_rich_vals[max(np.where(m_bin_count[:, i] > 0)[0])], 1.0, '*', c='black', markersize=10)
      ax1.set_xlim(min_rich_bin, max(m_rich_u_3))
      ax1.set_xscale("log", nonposx='clip')
      ax1.set_ylim(0, 1.4)
      ax1.set_xlabel('N$_{obs}$')
      ax1.set_ylabel('Completeness')
      string1 = 'Number of Haloes: ' + str(sum(m_bin_count[:, i]))
      string2 = 'Number of 1-Way Matches: ' + str(sum(m_match_bin_count_1[:, i]))
      string3 = 'Number of 2-Way Matches: ' + str(sum(m_match_bin_count_2[:, i]))
      ax1.text(0.05, 0.95, string1, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, size = 10)
      ax1.text(0.05, 0.9, string2, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, size = 10)
      ax1.text(0.05, 0.85, string3, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, size = 10)
      string4 = str(y_z_vals[i]) + ' > z > ' + str(y_z_vals[i + 1])  
      ax1.set_title(string4)
      ax1.legend()

      fig1_name = cluster_file + '_completeness_' + str(i).zfill(2) + '.pdf'
      fig1.savefig(fig1_name)
      print "Plots saved to:", fig1_name

      fig2 = plt.figure() 
      plt.rc('legend',**{'fontsize':6})
      gs2 = gridspec.GridSpec(1, 1)
      ax2 = fig2.add_subplot(gs2[0])

      majorFormatter = FormatStrFormatter('%4.2f')
      ax2.yaxis.set_major_formatter(majorFormatter)
      ax2.plot(xyline, (1, 1), '--', c='0.5', linewidth=1.0)
      ax2.plot(x_c_rich_vals, pure_bin_1[:, i], linewidth=1.5, c='b', label='One-Way')                #one-way purity plot
      ax2.plot(x_c_rich_vals, pure_bin_2[:, i], linewidth=1.5, c='r', label='Two-Way')                #two-way purity plot
      if (np.where(c_bin_count[:, i] > 0)[0]).size > 0:
         ax2.plot(x_c_rich_vals[max(np.where(c_bin_count[:, i] > 0)[0])], 1.0, '*', c='black', markersize=10)
      ax2.set_xlim(min_rich_bin, max(c_rich_u))
      ax2.set_xscale("log", nonposx='clip')
      ax2.set_ylim(0, 1.4)
      ax2.set_xlabel('N$_{obs}$')
      ax2.set_ylabel('Purity')
      string1 = 'Number of Clusters: ' + str(sum(c_bin_count[:, i]))
      string2 = 'Number of 1-Way Matches: ' + str(sum(c_match_bin_count_1[:, i]))
      string3 = 'Number of 2-Way Matches: ' + str(sum(c_match_bin_count_2[:, i]))
      ax2.text(0.05, 0.85, string1, horizontalalignment = 'left', verticalalignment = 'center', transform=ax1.transAxes, size = 10)
      ax2.text(0.05, 0.8, string2, horizontalalignment = 'left', verticalalignment = 'center', transform=ax1.transAxes, size = 10)
      ax2.text(0.05, 0.75, string3, horizontalalignment = 'left', verticalalignment = 'center', transform=ax1.transAxes, size = 10)
      string4 = str(y_z_vals[i]) + ' > z > ' + str(y_z_vals[i + 1])  
      ax2.set_title(string4)
      ax2.legend()

      fig2_name = cluster_file + '_purity_' + str(i).zfill(2) + '.pdf'
      fig2.savefig(fig2_name)
      print "Plots saved to:", fig2_name
   
   fig3 = plt.figure() 
   plt.hist(c_z_u, n_z_bins, range = [min_z_bin, max_z_bin], histtype = 'step', color = 'blue', label = 'Clusters')  #N(z) histogram plot
   plt.hist(m_z_u_3, n_z_bins, range = [min_z_bin, max_z_bin], histtype = 'step', color = 'red', linestyle = 'dashed', label = 'Haloes')
   plt.legend()
   plt.xlabel('$z$')
   plt.title('N(z) Distribution')

   fig3_name = cluster_file + '_nz.pdf'
   fig3.savefig(fig3_name)
   print "Plots saved to:", fig3_name
   
   print '================================================================='
   
if __name__ == "__main__":
   main(sys.argv[1:])
