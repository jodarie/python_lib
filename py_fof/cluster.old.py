###################
# CLUSTER CLASSES #
###################

import math, numpy as np
from functions import astro, comp, cosmo

class Cluster:
    """
    Class for storing cluster members and deriving properties from those members.
    """
    def __init__(self, c_id):
        """
        Function that initialises a Cluster instance.
        """
        self.id = c_id
        self.g_num = []
        self.g_id = []
        self.g_ra = []
        self.g_dec = []
        self.g_z = []
        self.match_flags = []
        self.match_flag = False
        self.match_list = []
    def clean(self, c_id):
        """
        Function that resets the cluster ID and empties the match flags.
        """
        self.id = c_id
        self.match_flags = []
        self.match_flag = False
        self.match_list = []
    def extend(self, g_num, g_id, g_ra, g_dec, g_z):
        """
        Function that adds galaxy properties to Cluster.
        """
        self.g_num.extend(g_num)
        self.g_id.extend(g_id)
        self.g_ra.extend(g_ra)
        self.g_dec.extend(g_dec)
        self.g_z.extend(g_z)
    def flag(self, c_id):
        """
        Function that adds Cluster ID to list of match flags.
        """
        self.match_flags.extend(c_id)
    def props(self, bg_expect):
        """
        Function that sets Cluster properties.
        """                
        self.ra = np.median(self.g_ra)
        self.dec = np.median(self.g_dec)
        self.z = np.median(self.g_z)
        self.ngal = len(self.g_id)
        distances = []        
        for i in range(self.ngal):
            distances.extend([astro.projected_distance(self.g_ra[i], self.ra,
                                                                 self.g_dec[i], self.dec)])
        self.size = np.mean(distances)
        self.area = self.size ** 2 * math.pi
        if bg_expect != None:
            self.sn = (self.ngal) / ((self.area * bg_expect) ** 0.5)
        else:
            self.sn = 0.0
    def unique(self):
        """
        Function that removes duplicate galaxies from Cluster.
        """
        unique_ids, index = np.unique(self.g_id, return_index = True)
        self.g_id = list(unique_ids)
        self.g_num = list(np.array(self.g_num)[index])
        self.g_ra = list(np.array(self.g_ra)[index])
        self.g_dec = list(np.array(self.g_dec)[index])
        self.g_z = list(np.array(self.g_z)[index])

class Galaxy:
    """
    Class for storing galaxy properties.
    """
    def __init__(self, id, ra, dec, z, dz, number):
        """
        Function that initialises a Galaxy instance.
        """
        self.num = number
        self.id = int(id)
        self.ra = float(ra)
        self.dec = float(dec)
        self.z = float(z)
        self.dz = float(dz)
        self.v = self.z / (1 + self.z)
    def assign_dis(self, c, H0, Omega_M, Omega_L):
        """
        Function that calculates the angular diameter distance
        for a Galaxy instance.
        """
        self.da = ((c / H0) * cosmo.angdidis(self.z, Omega_M, Omega_L))
    def assign_bin(self, bins, z_min, z_bin_size):
        """
        Function that assigns a Galaxy instance to the
        corresponding Z_bin instances.
        """
        self.bin = comp.find_bin(self.z, z_min, z_bin_size)        
        bins[self.bin].count += 1
    def assign_clt(self, nbins):
        self.clt = [None] * nbins        

class Z_bin:
    """
    Class for storing redshift bin properties.
    """
    def __init__(self, id, z, link_z):
        """
        Function that initialises a Z_bin instance.
        """
        self.id = id
        self.z = z
        self.count = 0
        self.r_ref = 0
        self.link_z = link_z / (1 + self.z)
    def assign_dis(self, c, H0, Omega_M, Omega_L):
        """
        Function that calculates the angular diameter distance
        for a Z_bin instance.
        """
        self.dcmv = cosmo.dcomvoldz(self.z, Omega_M, Omega_L)
        self.da = ((c / H0) * cosmo.angdidis(self.z, Omega_M, Omega_L))
    def linking_param(self, r_ref, z_bin_size):
        """
        Function that calculates the linking length for each
        redshift bin.
        """
        self.dndz = self.count / z_bin_size
        self.dndv = self.dndz / self.dcmv
        if self.dndv > 0:     
            self.link_r = self.dndv ** -0.5 * r_ref
        else:
            self.link_r = 0
        self.rfriend = self.link_r / self.da
