import numpy as np
from functions import astro

class KDTreeNode():
    def __init__(self, mem):
        self.mem = np.array(mem)
        self.size = len(self.mem)
        self.ra = np.median(self.mem[:, 0])
        self.dec = np.median(self.mem[:, 1])

        min_ra = min(self.mem[:, 0])
        min_dec = max(self.mem[:, 1])
        
        self.radius = astro.deg2rad(astro.projected_distance(self.ra, min_ra, self.dec, min_dec) / 60.0)
            
class KDTree():    
    def __init__(self, data, max_depth):

        self.node_list = []
        self.max_depth = max_depth
        
        def build_kdtree(point_list, depth):            
            axis = depth % 2
            point_list.sort(key = lambda point: point[axis])
            median = len(point_list) / 2
            if depth <= self.max_depth:
                if depth == self.max_depth:
                    node = KDTreeNode(point_list)
                    self.node_list.append(node)
                build_kdtree(point_list[0:median], depth + 1)
                build_kdtree(point_list[median:], depth + 1)
                
        build_kdtree(data, depth = 0)

        
