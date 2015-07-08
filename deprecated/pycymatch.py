#######################
#   PYCYMATCH V.3.2   #
#######################
# Samuel Farrens 2014 #
#######################

import sys, getopt, math, pycymatch_help as pch
import astro, comp, errors
import numpy as np, matplotlib as mp, scipy.stats as ss
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

def main(argv):
   
    ################
    # SET DEFAULTS #
    ################

    cluster_file = ''
    mock_file = ''
    delta_z = 0.03
    m_ra_col = 1 - 1
    m_dec_col = 2 - 1
    m_z_col = 3 - 1
    m_id_col = 4 - 1
    m_rich_col = 5 - 1
    m_mass_col = 6 - 1
    m_size_col = 8 - 1
    m_minra_col = 9 - 1
    m_maxra_col = 10 - 1
    m_mindec_col = 11 - 1
    m_maxdec_col = 12 - 1
    c_id_col = 1 - 1
    c_ra_col = 2 - 1
    c_dec_col = 3 - 1
    c_z_col = 4 - 1
    c_rich_col = 5 - 1
    c_sn_col = 6 - 1
    min_rich_bin = 0
    max_rich_bin = -1
    rich_bin_size = 0.2
    min_mass_bin = 13.0
    max_mass_bin = -1
    mass_bin_size = 0.1
    min_z_bin = 0.0
    max_z_bin = -1
    z_bin_size = 0.1
    
   ##################
   # READ ARGUMENTS #
   ##################
   
   try:
      opts, args = getopt.getopt(argv,"c:m:z:l:n:o:p:q:r:s:t:u:v:w:x:a:b:d:e:f:g:i:j:k:h",["cfile=","mfile=","delta_z=","m_id_col=","m_ra_col=",
                                                                                         "m_dec_col=","m_z_col=","m_rich_col=","m_mass_col=",
                                                                                         "m_size_col=","c_ra_col=","c_dec_col=","c_z_col=",
                                                                                         "c_rich_col=","c_sn_col=","min_rich_bin=","max_rich_bin=",
                                                                                         "rich_bin_size=","min_mass_bin=","max_mass_bin=",
                                                                                         "mass_bin_size=","min_z_bin=","max_z_bin=","z_bin_size=",
                                                                                         "help"])
   except getopt.GetoptError:
        pch.print_help()
   for opt, arg in opts:
      if opt in ('-h','--help'):
            pch.print_help()
      elif opt in ("-c", "--cfile"):                                                       #set cluster file name
         cluster_file = arg
      elif opt in ("-m", "--mfile"):                                                       #set mock halo file name
         mock_file = arg
      elif opt in ("-z", "--delta_z"):                                                     #set redshift length of cylinder
         delta_z = float(arg)
      elif opt in ("-l", "--m_id_col"):                                                    #set mock halo id column number
         m_id_col = int(arg) - 1
      elif opt in ("-n", "--m_ra_col"):                                                    #set mock halo ra column number
         m_ra_col = int(arg) - 1
      elif opt in ("-o", "--m_dec_col"):                                                   #set mock halo dec column number
         m_dec_col = int(arg) - 1
      elif opt in ("-p", "--m_z_col"):                                                     #set mock halo z column number
         m_z_col = int(arg) - 1
      elif opt in ("-q", "--m_rich_col"):                                                  #set mock halo richness column number
         m_rich_col = int(arg) - 1
      elif opt in ("-r", "--m_mass_col"):                                                  #set mock halo mass column number
         m_mass_col = int(arg) - 1
      elif opt in ("-s", "--m_size_col"):                                                  #set mock halo size column number
         m_size_col = int(arg) - 1
      elif opt in ("-t", "--c_ra_col"):                                                    #set cluster ra column number
         c_ra_col = int(arg) - 1
      elif opt in ("-u", "--c_dec_col"):                                                   #set cluster dec column number
         c_dec_col = int(arg) - 1
      elif opt in ("-v", "--c_z_col"):                                                     #set cluster z column number
         c_z_col = int(arg) - 1
      elif opt in ("-w", "--c_rich_col"):                                                  #set cluster rich column number
         c_rich_col = int(arg) - 1
      elif opt in ("-x", "--c_sn_col"):                                                    #set cluster S/N column number
         c_sn_col = int(arg) - 1
      elif opt in ("-a", "--min_rich_bin"):                                                #set minimum for binning in richness
         min_rich_bin = float(arg)
      elif opt in ("-b", "--max_rich_bin"):                                                #set maximum for binning in richness
         max_rich_bin = float(arg)
      elif opt in ("-d", "--rich_bin_size"):                                               #set size of bins in richness
         rich_bin_size = float(arg)
      elif opt in ("-e", "--min_mass_bin"):                                                #set minimum for binning in mass
         min_mass_bin = float(arg)
      elif opt in ("-f", "--max_mass_bin"):                                                #set maximum for binning in mass
         max_mass_bin = float(arg)
      elif opt in ("-g", "--mass_bin_size"):                                               #set size of bins in mass
         mass_bin_size = float(arg)
      elif opt in ("-i", "--min_z_bin"):                                                   #set minimum for binning in redshift
         min_z_bin = float(arg)
      elif opt in ("-j", "--max_z_bin"):                                                   #set maximum for binning in redshift
         max_z_bin = float(arg)
      elif opt in ("-k", "--z_bin_size"):                                                  #set size of bins in redshift
         z_bin_size = float(arg)
      else:
         print 'Something is wrong!'
         sys.exit(2)

   #############
   # READ DATA #
   #############
         
   errors.file_name_error(mock_file)                                                       #make sure files exit
   errors.file_name_error(cluster_file)

   print '================================================================='
   
   print 'Reading file:', mock_file,        
   mock_data = np.genfromtxt(mock_file, dtype="S", unpack = True)                          #mock haloes
   print '\tComplete'
   print 'Reading file:', cluster_file,        
   cluster_data = np.genfromtxt(cluster_file, dtype="S", unpack = True)                    #clusters
   print '\tComplete'

   ####################
   # DEFINE VARIABLES #
   ####################

   m_id = np.array(mock_data[m_id_col, :])
   m_ra = np.array(mock_data[m_ra_col, :], dtype = 'f')
   m_dec = np.array(mock_data[m_dec_col, :], dtype = 'f')
   m_z = np.array(mock_data[m_z_col, :], dtype = 'f')
   m_mass = np.array(mock_data[m_mass_col, :], dtype = 'f')
   m_rich = np.log10(np.array(mock_data[m_rich_col, :], dtype = 'f'))
   m_size = np.array(mock_data[m_size_col, :], dtype = 'f')
   c_id = np.array(cluster_data[c_id_col, :], dtype = 'i')
   c_ra = np.array(cluster_data[c_ra_col, :], dtype = 'f')
   c_dec = np.array(cluster_data[c_dec_col, :], dtype = 'f')
   c_z = np.array(cluster_data[c_z_col, :], dtype = 'f')
   c_rich = np.log10(np.array(cluster_data[c_rich_col, :], dtype = 'f'))
   c_sn = np.array(cluster_data[c_sn_col, :], dtype = 'f')

   #######################
   # SET BOUNDARY LIMITS #
   #######################

   if max_rich_bin < 0: max_rich_bin = max(max(c_rich), max(m_rich))
   if max_mass_bin < 0: max_mass_bin = max(m_mass)
   if max_z_bin < 0: max_z_bin = max(max(c_z), max(m_z))

   m_limit_index = ((m_z >= min_z_bin) & (m_z <= max_z_bin) &
                    (m_rich >= min_rich_bin) & (m_rich <= max_rich_bin) &
                    (m_mass >= min_mass_bin) & (m_mass <= max_mass_bin))
   c_limit_index = ((c_z >= min_z_bin) & (c_z <= max_z_bin) &
                    (c_rich >= min_rich_bin) & (c_rich <= max_rich_bin))

   m_id = m_id[m_limit_index]
   m_ra = m_ra[m_limit_index]
   m_dec = m_dec[m_limit_index]
   m_z = m_z[m_limit_index]
   m_mass = m_mass[m_limit_index]
   m_rich = m_rich[m_limit_index]
   m_size = m_size[m_limit_index]
   c_id = c_id[c_limit_index]
   c_ra = c_ra[c_limit_index]
   c_dec = c_dec[c_limit_index]
   c_z = c_z[c_limit_index]
   c_rich = c_rich[c_limit_index]
   c_sn = c_sn[c_limit_index]

   ####################
   # SORT BY RICHNESS #
   ####################
      
   m_rich_index = m_rich.argsort()[::-1]
   c_rich_index = c_rich.argsort()[::-1]

   m_id = m_id[m_rich_index]
   m_ra = m_ra[m_rich_index]
   m_dec = m_dec[m_rich_index]
   m_z = m_z[m_rich_index]
   m_mass = m_mass[m_rich_index]
   m_rich = m_rich[m_rich_index]
   m_size = m_size[m_rich_index]
   c_id = c_id[c_rich_index]
   c_ra = c_ra[c_rich_index]
   c_dec = c_dec[c_rich_index]
   c_z = c_z[c_rich_index]
   c_rich = c_rich[c_rich_index]
   c_sn = c_sn[c_rich_index]
   
   ###################
   # BIN BY RICHNESS #
   ###################

   n_rich_bins = int(math.floor((max_rich_bin - min_rich_bin) / rich_bin_size)) + 1 #set number of richness bins

   c_rich_bin_index = np.floor((c_rich - min_rich_bin) / rich_bin_size).astype('int')
   m_rich_bin_index = np.floor((m_rich - min_rich_bin) / rich_bin_size).astype('int')

   x_rich_vals = (np.arange(n_rich_bins) + 0.5) * rich_bin_size + min_rich_bin

   ############
   # BIN BY Z #
   ############

   n_z_bins = int(math.floor((max_z_bin - min_z_bin) / z_bin_size)) + 1 #set number of redshift bins

   c_z_bin_index = np.floor((c_z - min_z_bin) / z_bin_size).astype('int')
   m_z_bin_index = np.floor((m_z - min_z_bin) / z_bin_size).astype('int')

   x_z_vals = (np.arange(n_z_bins) + 0.5) * z_bin_size + min_z_bin

   ################
   # FIND MATCHES #
   ################

   delta_angle = 1.0 #to speed up calculations
   z_factor = 2.0 #to adjust the line-of-sight matching
   r_factor = 0.5 #to adjust the transverse matching 
   
   c_bin_count = np.zeros((n_rich_bins, n_z_bins)).astype('int')
   m_bin_count = np.zeros((n_rich_bins, n_z_bins)).astype('int')
   c_match_bin_count = np.zeros((n_rich_bins, n_z_bins)).astype('int')
   m_match_bin_count = np.zeros((n_rich_bins, n_z_bins)).astype('int')

   m_match_flag = np.zeros(len(m_ra)).astype('int') - 1
   c_match_flag = np.zeros(len(c_ra)).astype('int') - 1

   for i in range(len(c_ra)):                                                          
      c_bin_count[c_rich_bin_index[i], c_z_bin_index[i]] += 1
   for i in range(len(m_ra)):
      m_bin_count[m_rich_bin_index[i], m_z_bin_index[i]] += 1
      index = np.where((c_match_flag == -1) & (np.fabs(m_z[i] - c_z) <= (z_factor * delta_z * (1 + m_z[i]))) &
                      (np.fabs(m_ra[i] - c_ra) <= delta_angle) & (np.fabs(m_dec[i] - c_dec) <= delta_angle))[0]
      for j in index:
         if m_match_flag[i] == -1:
            dist = astro.projected_distance(m_ra[i], c_ra[j], m_dec[i], c_dec[j])
            if dist <= (m_size[i] * r_factor):
               c_match_bin_count[c_rich_bin_index[j], c_z_bin_index[j]] += 1
               m_match_bin_count[m_rich_bin_index[i], m_z_bin_index[i]] += 1               
               m_match_flag[i] = j
               c_match_flag[j] = i
               
   pure_bin = np.flipud(np.array(c_match_bin_count).astype("float") / np.array(c_bin_count).astype("float"))
   complete_bin = np.flipud(np.array(m_match_bin_count).astype("float") / np.array(m_bin_count).astype("float"))
     
   match_index = np.where(m_match_flag > -1)[0]

   #############################
   # PRINT M > m_limit MATCHES #
   #############################

   m_limit = 13.0 #mass limit
   
   match_out_file = cluster_file + '_halo_matching.txt'                                          
   match_out = open(match_out_file,'w')

   mass_sort_index = m_mass.argsort()[m_mass[m_mass.argsort()] >= m_limit][::-1]   
      
   print>> match_out, '#H_ID                  H_RA   H_DEC  H_Z   H_MASS H_NGAL MATCH   C_RA   C_DEC  C_Z   C_NGAL C_SN'

   for i in mass_sort_index:
    print>> match_out, '%-22s' % m_id[i],'%06.3f' % m_ra[i], '%+06.3f' % m_dec[i], '%05.3f' % m_z[i],
    print>> match_out, '%06.3f' % m_mass[i], '%04i' % (10 ** m_rich[i]),
    if m_match_flag[i] > -1:
        print>> match_out, '  Y      ', '%06.3f' % c_ra[m_match_flag[i]], '%+06.3f' % c_dec[m_match_flag[i]],
        print>> match_out, '%05.3f' % c_z[m_match_flag[i]], '%04i  ' % (10 ** c_rich[m_match_flag[i]]),
        print>> match_out, '%06.3f' % c_sn[m_match_flag[i]] 
    else:
        print>> match_out, '  N'

   print "Matching data saved to:", match_out_file

   #########################
   # PRINT CLUSTER MATCHES #
   #########################

   c_match_out_file = cluster_file + '_cluster_matching.txt'                                          
   c_match_out = open(c_match_out_file,'w')
      
   print>> c_match_out, '#C_RA  C_DEC  C_Z   C_NGAL C_S/N    MATCH   H_ID                   H_RA   H_DEC  H_Z   H_NGAL H_MASS'

   for i in range(len(c_rich)):
    print>> c_match_out, '%06.3f' % c_ra[i], '%+06.3f' % c_dec[i], '%05.3f' % c_z[i],
    print>> c_match_out, '%04i' % (10 ** c_rich[i]), '  %06.3f' % c_sn[i],
    if c_match_flag[i] > -1:
        print>> c_match_out, '  Y      ', '%-22s' % m_id[c_match_flag[i]], '%06.3f' % m_ra[c_match_flag[i]], '%+06.3f' % m_dec[c_match_flag[i]],
        print>> c_match_out, '%05.3f' % m_z[c_match_flag[i]], '%04i  ' % (10 ** m_rich[c_match_flag[i]]),
        print>> c_match_out, '%06.3f' % m_mass[c_match_flag[i]] 
    else:
        print>> c_match_out, '  N'

   print "Matching data saved to:", c_match_out_file
  
   ###############
   # BIN BY MASS #
   ###############

   n_mass_bins = int(math.floor((max_mass_bin - min_mass_bin) / mass_bin_size)) + 1

   m_mass_bin_index = np.floor((m_mass - min_mass_bin) / mass_bin_size).astype('int')

   y_mass_vals = (np.arange(n_mass_bins) + 0.5) * mass_bin_size + min_mass_bin

   m_bin_count = np.zeros((n_mass_bins, n_z_bins)).astype('int')
   m_match_bin_count = np.zeros((n_mass_bins, n_z_bins)).astype('int')

   for i in range(len(m_mass)):
      m_bin_count[m_mass_bin_index[i], m_z_bin_index[i]] += 1
      if m_match_flag[i] > -1:
        m_match_bin_count[m_mass_bin_index[i], m_z_bin_index[i]] += 1               
               
   complete_bin_mass = np.flipud(np.array(m_match_bin_count).astype("float") / np.array(m_bin_count).astype("float"))
     
   ################
   # SPEARMAN RHO #
   ################
   
   for i in range(n_z_bins):
    m_index_1 = m_rich[match_index][m_z_bin_index[match_index] == i]
    m_index_2 = m_rich[match_index][c_z_bin_index[m_match_flag[match_index]] == i]
    c_index_1 = c_rich[m_match_flag[match_index]][m_z_bin_index[match_index] == i]
    c_index_2 = c_rich[m_match_flag[match_index]][c_z_bin_index[m_match_flag[match_index]] == i]
    if m_index_1.size > 1:
        m_rho, p1 = ss.spearmanr(m_index_1, c_index_1)
        m_rho_err = 0.6325 / (len(m_index_1) - 1) ** 0.5
    else:
        m_rho = 0.0
        m_rho_err = 0.0
    if m_index_2.size > 1:
        c_rho, p2 = ss.spearmanr(m_index_2, c_index_2)
        c_rho_err = 0.6325 / (len(m_index_2) - 1) ** 0.5
    else:
        c_rho = 0.0
        c_rho_err = 0.0
    if i == 0:
        m_rhos = np.array(m_rho)
        m_rhos_err = np.array(m_rho_err)
        c_rhos = np.array(c_rho)
        c_rhos_err = np.array(c_rho_err)
    else:
        m_rhos = np.hstack((m_rhos, m_rho))
        m_rhos_err = np.hstack((m_rhos_err, m_rho_err))
        c_rhos = np.hstack((c_rhos, c_rho))
        c_rhos_err = np.hstack((c_rhos_err, c_rho_err))
        
   ##############
   # MAKE PLOTS #
   ##############

   interpolation = 'spline16'
   colormap = cm.jet
   boundaries = np.arange(0.0, 1.1, 0.1)
   normalisation = colors.BoundaryNorm(boundaries, colormap.N)
   
   #Completeness vs. z and N_true_mock
   fig1 = plt.figure()
   fig1.subplots_adjust(wspace = .01)
   gs1 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
   ax11 = fig1.add_subplot(gs1[0])
   cax11 = fig1.add_subplot(gs1[1])
   im = ax11.imshow(complete_bin, extent = [min_z_bin, max_z_bin, min_rich_bin, max_rich_bin], aspect = 'auto',
                     interpolation = interpolation, cmap = colormap, norm = normalisation)
   ax11.set_xlabel('z')
   ax11.set_ylabel('log$_{10}$ N$_{true}$')
   cbar = fig1.colorbar(im, cax=cax11)
   cbar.set_label('Completeness')
   fig1_name = cluster_file + '_completeness_plot.pdf'
   fig1.savefig(fig1_name)
   print "Plots saved to:", fig1_name

   #Purity vs. z and N_obs
   fig2 = plt.figure()
   fig2.subplots_adjust(wspace = .01)
   gs2 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
   ax21 = fig2.add_subplot(gs2[0])
   cax21 = fig2.add_subplot(gs2[1])
   im = ax21.imshow(pure_bin, extent = [min_z_bin, max_z_bin, min_rich_bin, max_rich_bin], aspect = 'auto',
                     interpolation = interpolation, cmap = colormap, norm = normalisation)
   ax21.set_xlabel('z')
   ax21.set_ylabel('log$_{10}$ N$_{obs}$')
   cbar = fig2.colorbar(im, cax=cax21)
   cbar.set_label('Purity')
   fig2_name = cluster_file + '_purity_plot.pdf'
   fig2.savefig(fig2_name)
   print "Plots saved to:", fig2_name

   #Mass vs. N_true_mock
   heat, yedges, xedges = np.histogram2d(m_rich, m_mass, bins=(n_rich_bins, n_mass_bins))
   heat = np.log10(heat)   
   fig3 = plt.figure()
   fig3.subplots_adjust(wspace = .01)
   gs3 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
   ax31 = fig3.add_subplot(gs3[0])
   cax31 = fig3.add_subplot(gs3[1])
   im = ax31.imshow(np.flipud(heat), extent = [min_mass_bin, max_mass_bin, min_rich_bin, max_rich_bin], aspect = 'auto',
                     interpolation = interpolation, cmap = colormap)
   ax31.set_xlabel('log$_{10}$ M')
   ax31.set_ylabel('log$_{10}$ N$_{true}$')
   cbar = fig3.colorbar(im, cax=cax31)
   cbar.set_label('log$_{10}$ N')
   fig3_name = cluster_file + '_mock_plot1.pdf'
   fig3.savefig(fig3_name)
   print "Plots saved to:", fig3_name

   #z vs. N_true_mock
   heat, yedges, xedges = np.histogram2d(m_rich, m_z, bins=(n_rich_bins, n_z_bins))
   heat = np.log10(heat)   
   fig4 = plt.figure()
   fig4.subplots_adjust(wspace = .01)
   gs4 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
   ax41 = fig4.add_subplot(gs4[0])
   cax41 = fig4.add_subplot(gs4[1])
   im = ax41.imshow(np.flipud(heat), extent = [min_z_bin, max_z_bin, min_rich_bin, max_rich_bin], aspect = 'auto',
                     interpolation = interpolation, cmap = colormap)
   ax41.set_xlabel('z')
   ax41.set_ylabel('log$_{10}$ N$_{true}$')
   cbar = fig4.colorbar(im, cax=cax41)
   cbar.set_label('log$_{10}$ N')
   fig4_name = cluster_file + '_mock_plot2.pdf'
   fig4.savefig(fig4_name)
   print "Plots saved to:", fig4_name
   
   #z vs. spearman rho   
   fig5 = plt.figure()
   gs5 = gridspec.GridSpec(1, 1)
   ax51 = fig5.add_subplot(gs5[0])
   ax51.plot(x_z_vals, m_rhos, linewidth=1.0, c='b', label='z$_{mock}$ Bins')
   ax51.plot(x_z_vals, c_rhos, linewidth=1.0, linestyle='dashed', c='r', label='z$_{obs}$ Bins')
   ax51.errorbar(x_z_vals, m_rhos, yerr=m_rhos_err, fmt='bx')
   ax51.errorbar(x_z_vals, c_rhos, yerr=c_rhos_err, fmt='rx')
   ax51.set_ylim(0.0, 1.0)
   ax51.set_xlabel('z')
   ax51.set_ylabel(r'$\rho$')
   ax51.set_title('Spearman Rank Order Coeffcient')
   ax51.legend()
   fig5_name = cluster_file + '_spearman.pdf'
   fig5.savefig(fig5_name)
   print "Plots saved to:", fig5_name
   
   #Completeness vs. z and Mass
   fig6 = plt.figure()
   fig6.subplots_adjust(wspace = .01)
   gs6 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
   ax61 = fig6.add_subplot(gs6[0])
   cax61 = fig6.add_subplot(gs6[1])
   im = ax61.imshow(complete_bin_mass, extent = [min_z_bin, max_z_bin, min_mass_bin, max_mass_bin], aspect = 'auto',
                     interpolation=interpolation, cmap = colormap, norm = normalisation)
   ax61.set_xlabel('z')
   ax61.set_ylabel('log$_{10}$ M')
   cbar = fig6.colorbar(im, cax=cax61)
   cbar.set_label('Completeness')
   fig6_name = cluster_file + '_completeness_mass_plot.pdf'
   fig6.savefig(fig6_name)
   print "Plots saved to:", fig6_name

  #########################################################################
   
if __name__ == "__main__":
   main(sys.argv[1:])
