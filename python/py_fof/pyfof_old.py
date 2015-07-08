#######################
#     PYFOF V.1.0     #
#######################
# Samuel Farrens 2014 #
#######################

import math, optparse
from functions import astro, comp, cosmo, errors, merge, options
import numpy as np, scipy
from scipy import spatial
from classes.cluster import Cluster

# Classes:

class Galaxy:
    """
    Class for storing galaxy properties.
    """
    def __init__(self, data, number):
        """
        Function that initialises a Galaxy instance.
        """
        self.num = number
        self.id = data[opts.id_col - 1]
        self.ra = float(data[opts.ra_col - 1])
        self.dec = float(data[opts.dec_col - 1])
        self.z = float(data[opts.z_col - 1])
        self.dz = float(data[opts.dz_col - 1])
        self.v = self.z / (1 + self.z)
        self.da = ((opts.c / opts.H0) * cosmo.angdidis(self.z, opts.Omega_M, opts.Omega_L))
        self.x = math.sin(astro.deg2rad(90 - self.dec)) * math.cos(astro.deg2rad(self.ra)) * self.da
        self.y = math.sin(astro.deg2rad(90 - self.dec)) * math.sin(astro.deg2rad(self.ra)) * self.da
        self.photoz_bins = []
    def assign_bin(self, bins):
        """
        Function that assigns a Galaxy instance to the
        corresponding Z_bin instances.
        """
        self.bin = comp.find_bin(self.z, opts.z_min, opts.z_bin_size)        
        bins[self.bin].count += 1
        if opts.mode == 'phot':
            min_z = comp.find_bin(round(self.z - self.dz * opts.link_z, 2), opts.z_min, opts.z_bin_size) #can be changed
            max_z = comp.find_bin((self.z + self.dz * opts.link_z), opts.z_min, opts.z_bin_size)            
            num_bins = comp.num_bins(min_z, max_z, 1) + 1
            for i in range(num_bins):
                self.photoz_bins.append(min_z + i)
    def assign_clt(self, nbins):
        self.clt = [None] * nbins        

class Z_bin:
    """
    Class for...
    """
    def __init__(self, id, z):
        """
        Function that initialises a Z_bin instance.
        """
        self.id = id
        self.z = z
        self.count = 0
        self.r_ref = 0
        self.dcmv = cosmo.dcomvoldz(z, opts.Omega_M, opts.Omega_L)
        self.link_z = opts.link_z / (1 + self.z)
        self.da = ((opts.c / opts.H0) * cosmo.angdidis(self.z, opts.Omega_M, opts.Omega_L))
    def linking_param(self, r_ref):
        """
        Function that calculates the linking length for each
        redshift bin.
        """
        self.dndz = self.count / opts.z_bin_size
        self.dndv = self.dndz / self.dcmv
        if self.dndv > 0:     
            self.link_r = self.dndv ** -0.5 * r_ref
        else:
            self.link_r = 0
        self.rfriend = self.link_r / self.da
            
# Functions:

def set_up_z_bins(z_min, z_max, z_bin_size):
    """
    Function that sets up the redshift bins.
    """
    n_z_bins = int(((z_max - z_min) / z_bin_size) + 1)
    z_bins = []
    for i in range(n_z_bins):
        z_bins.append(Z_bin(i, z_min + i * z_bin_size))    
    return np.array(z_bins)

def assign_linking_parameters(z_bins):
    """
    Function that finds which redshift bin each galaxy
    corresponds to.
    """
    z_ref_index = np.floor(((opts.z_ref - opts.z_min) / opts.z_bin_size) + 0.5).astype('int')
    r_ref = (z_bins[z_ref_index].count / (opts.z_bin_size * z_bins[z_ref_index].dcmv)) ** 0.5 * opts.link_r
    for i in range(z_bins.size):
        z_bins[i].linking_param(r_ref)

def friendship(zbin, gal1, gal2, mode):
    """
    Function that checks if two galaxies are friends in
    a given redshift bin.
    """
    if mode == 'spec':
        rfriend = zbin.link_r / gal1.da
    else:
        rfriend = zbin.rfriend
    dist =  astro.deg2rad(astro.projected_distance(gal1.ra, gal2.ra, gal1.dec, gal2.dec) / 60.0)
    check1 = (dist < rfriend)

    if gal1.id == '7418' and gal2.id == '9651':
        print dist, rfriend, check1
        if check1:
            print zbin.z, gal2.id
    
    if mode == 'spec':
        check2 = (math.fabs(gal1.v - gal2.v) <= zbin.link_z)
        check = (check1 & check2)
    else:
        check = check1
    return check

def bin_check(mode, gal_bins, bin):
    """
    Function that checks if a bla ....
    """
    if mode == 'spec':
        check = True
    else:
        check = bin in gal_bins
    return check

def bin_check2(mode, gal, zbin):
    """
    Function that checks if a bla ....
    """
    if mode == 'spec':
        check = True
    else:
        if np.abs(gal.z - zbin.z) <= opts.link_z * gal.dz:
            check = True
        else:
            check = False
    return check

def friends_of_friends(zbin, gals, kdtree, tree_dist, mode, cluster_count):
    """
    Function that looks for galaxy friends in a given
    redshift bin and returns a list of the resulting
    Cluster instances.
    """
    #loop over galaxies
    cluster_list = []
    for i in range(len(gals)):
        if  mode == 'spec':
            link_zbin = z_bins[gals[i].bin]
        else:
            link_zbin = z_bins[zbin]
        if bin_check2(mode, gal[i], link_zbin):
        #if bin_check(mode, gal[i].photoz_bins, zbin):
            #friends loop
            if gals[i].clt[zbin] == None:
                list = kdtree.query_ball_point([gals[i].x, gals[i].y], tree_dist)
                for j in list:
                    if ((bin_check2(mode, gal[j], link_zbin)) & (gals[i].id != gals[j].id) & (gals[j].clt[zbin] == None)):
                    #if ((bin_check(mode, gal[j].photoz_bins, zbin)) & (gals[i].id != gals[j].id) & (gals[j].clt[zbin] == None)):
                        if friendship(link_zbin, gals[i], gals[j], opts.mode):

                            if zbin == 26 and gal[i].id == '2353':
                                print gal[j].id

                              
                            if gals[i].clt[zbin] == None:
                                gals[i].clt[zbin] = len(cluster_list)
                                gals[j].clt[zbin] = gals[i].clt[zbin]
                                cluster_list.append(Cluster(cluster_count))
                                cluster_count += len(cluster_list)
                                cluster_list[gals[i].clt[zbin]].extend([gals[i].num], [gals[i].id], [gals[i].ra],
                                                                        [gals[i].dec], [gals[i].z], [gals[i].x], [gals[i].y])
                                cluster_list[gals[i].clt[zbin]].extend([gals[j].num], [gals[j].id], [gals[j].ra],
                                                                        [gals[j].dec], [gals[j].z], [gals[j].x], [gals[j].y])
                            else:
                                gals[j].clt[zbin] = gals[i].clt[zbin]
                                cluster_list[gals[i].clt[zbin]].extend([gals[j].num], [gals[j].id], [gals[j].ra],
                                                                        [gals[j].dec], [gals[j].z], [gals[j].x], [gals[j].y])
            #friends-of-friends loop
            if gals[i].clt[zbin] != None:
                clt_now = cluster_list[gals[i].clt[zbin]]
                for k in range(len(clt_now.g_id)):
                    if(gals[i].id != clt_now.g_id[k]):
                        list = kdtree.query_ball_point([clt_now.g_x[k], clt_now.g_y[k]], tree_dist)
                        for l in list:
                            if ((bin_check2(mode, gal[l], link_zbin)) & (clt_now.g_id[k] != gals[l].id)
                                & (gals[l].clt[zbin] == None)):
                            #if ((bin_check(mode, gal[l].photoz_bins, zbin)) & (clt_now.g_id[k] != gals[l].id)
                              #  & (gals[l].clt[zbin] == None)):
                                if friendship(link_zbin, gals[clt_now.g_num[k]], gals[l], opts.mode):
                                    gals[l].clt[zbin] = gals[i].clt[zbin]
                                    cluster_list[gals[i].clt[zbin]].extend([gals[l].num], [gals[l].id], [gals[l].ra],
                                                                        [gals[l].dec], [gals[l].z], [gals[l].x], [gals[l].y])
    return cluster_list

def assign_cluster_props(cluster_list):
    """
    Funcion that assigns Cluster properties and
    deletes those with less than the minimum number
    of members.
    """
    for i in range(len(cluster_list) - 1, -1, -1):
        cluster_list[i].unique()
        cluster_list[i].props(None)
        if cluster_list[i].ngal < opts.min_mem:
            del cluster_list[i]
    return cluster_list

##################
# READ ARGUMENTS #
##################

parser = optparse.OptionParser()
options.single_input(parser)
options.single_output(parser)
options.cosmo(parser)
options.fof(parser)
options.merge(parser)
(opts, args) = parser.parse_args()
  
if not opts.input_file:
    parser.error('Input file name not provided.')
if not opts.output_file:
    opts.output_file = opts.input_file + '_clusters.dat'
if not opts.link_r:
    parser.error('Transverse linking length not provided.')
if not opts.link_z:
    parser.error('Line-of-sight linking length not provided.')
    
#############
# READ DATA #
#############

#make sure files exit
errors.file_name_error(opts.input_file)

#read input file
print 'Reading file:', opts.input_file
data = np.genfromtxt(opts.input_file, unpack = True, dtype="S")

################
# SET DEFAULTS #
################

all_gal_z = np.array(data[opts.z_col - 1]).astype('float')

#set redshift limits
if not opts.z_min or not opts.z_max:
    if not opts.z_min: opts.z_min = np.min(all_gal_z)
    if not opts.z_max: opts.z_max = np.max(all_gal_z)

#restrict data to redshift limits
data = data[:, ((all_gal_z >= opts.z_min) & (all_gal_z < opts.z_max))]

del all_gal_z

#set up redshift bins
z_bins = set_up_z_bins(opts.z_min, opts.z_max, opts.z_bin_size)

############################
# DEFINE GALAXY PROPERTIES #
############################

#assign input data to Galaxy class
print 'Assign galaxy properties...'
gal = []
for i in range(len(data[0])):
    gal.append(Galaxy(data[:,i], i))
    gal[i].assign_bin(z_bins)
    if opts.mode == 'phot':
        gal[i].assign_clt(z_bins.size)
    else:
        gal[i].assign_clt(1)

################
# MAKE KD-TREE #
################

#make spatial kd-tree of galaxy positions
pos = []
for i in gal:
    pos.append([i.x, i.y])
tree = scipy.spatial.cKDTree(pos)
del pos

############################
# SET LININKING PARAMETERS #
############################

#set up bin properties and assign galaxies to bins
print 'Assign linking parameters...'
assign_linking_parameters(z_bins)

######################
# SEARCH FOR FRIENDS #
######################

if opts.mode == 'phot':
    loop_over_zbins = len(z_bins)
else:
    loop_over_zbins = 1

#loop through redshift bins and perform friends-of-friends
print 'Find clusters...'
clusters = []
for i in range(loop_over_zbins):
    clusters.extend(friends_of_friends(i, gal, tree, opts.tree_dist, opts.mode, len(clusters)))
    print '>>', z_bins[i].z, len(clusters)
    
#assign properties to detected clusters
assign_cluster_props(clusters)

print 'Found', len(clusters), 'clusters.'

##################
# MERGE CLUSTERS #
##################

#merge clusters with galaxies in common
print 'Merge clusters...'
merge.merge_clusters(clusters, opts.progress, opts.bg_expect, 1.0, 0.2)

print 'Merged into', len(clusters), 'clusters.'

###################
# OUTPUT CLUSTERS #
###################

ngal_list = []
for x in clusters:
    ngal_list.append(x.ngal)
index = np.array(ngal_list).argsort()[::-1]

clusters = np.array(clusters)[index]

output = open(opts.output_file, 'w')

for x in clusters:
    for y in range(x.ngal):
        print>> output, x.id, x.ngal, x.g_id[y], x.g_ra[y], x.g_dec[y], x.g_z[y]
