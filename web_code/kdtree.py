from functions import astro

""" Angular KD-Tree Class

Modified from code produced by
Matej Drame [matej.drame@gmail.com]
"""

class KDTreeNode():
    def __init__(self, point, left, right):
        self.point = point
        self.left = left
        self.right = right
    
    def is_leaf(self):
        return (self.left == None and self.right == None)

class KDTreeNeighbours():
    """ Internal structure used in nearest-neighbours search.
    """
    def __init__(self, query_point, radius):
        self.query_point = query_point
        self.radius = radius
        self.current_best = []

    def add(self, point):
        #angular separation in arcminutes
        ang_sep = astro.projected_distance(point[0], self.query_point[0], point[1], self.query_point[1])
        if ang_sep <= self.radius:
            self.current_best.append([point[2]])
    
    def get_best(self):
         return self.current_best
    #    return [element[0] for element in self.current_best]

        
class KDTree():
    """ KDTree implementation.
    
        Example usage:
        
            from kdtree import KDTree
            
            data = <load data> # iterable of points (which are also iterable, same length)
            point = <the point of which neighbours we're looking for>
            
            tree = KDTree.construct_from_data(data)
            nearest = tree.query(point, 30) # find points within 30 arcminutes
    """
    
    def __init__(self, data):
        def build_kdtree(point_list, depth):            
            # code based on wikipedia article: http://en.wikipedia.org/wiki/Kd-tree
            if not point_list:
                return None

            # select axis based on depth so that axis cycles through all valid values
            axis = depth % 2
            
            # sort point list and choose median as pivot point,
            point_list.sort(key=lambda point: point[axis])
            
            median = len(point_list)/2 # choose median

            # create node and recursively construct subtrees
            node = KDTreeNode(point=point_list[median],
                              left=build_kdtree(point_list[0:median], depth+1),
                              right=build_kdtree(point_list[median+1:], depth+1))
            return node
        
        self.root_node = build_kdtree(data, depth=0)
    
    @staticmethod
    def construct_from_data(data):
        tree = KDTree(data)
        return tree

    def query(self, query_point, radius):
        
        def nn_search(node, query_point, radius, depth, best_neighbours):
            if node == None:
                return
                        
            # if we have reached a leaf, let's add to current best neighbours,
            # (if it's better than the worst one or if there is not enough neighbours)
            if node.is_leaf():
                best_neighbours.add(node.point)
                return
            
            # this node is no leaf
            
            # select dimension for comparison (based on current depth)
            axis = depth % 2
            
            # figure out which subtree to search
            near_subtree = None # near subtree
            far_subtree = None # far subtree (perhaps we'll have to traverse it as well)
            
            # compare query_point and point of current node in selected dimension
            # and figure out which subtree is farther than the other
            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left

            # recursively search through the tree until a leaf is found
            nn_search(near_subtree, query_point, radius, depth+1, best_neighbours)

            # while unwinding the recursion, check if the current node
            # is closer to query point than the current best,
            # also, until t points have been found, search radius is infinity
            best_neighbours.add(node.point)
            
            # check whether there could be any points on the other side of the
            # splitting plane that are closer to the query point than the current best
            if astro.projected_distance(node.point[0], query_point[0], node.point[1], query_point[1]) <= radius:
                nn_search(far_subtree, query_point, radius, depth+1, best_neighbours)
    
            return
        
        # if there's no tree, there's no neighbors
        if self.root_node != None:
            neighbours = KDTreeNeighbours(query_point, radius)
            nn_search(self.root_node, query_point, radius, depth=0, best_neighbours=neighbours)
            result = neighbours.get_best()
        else:
            result = []
        
        return result
