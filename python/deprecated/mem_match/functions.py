import numpy as np

def find_bin(val, min_val, bin_size):
    '''
    Fine corresponding bin of a given value.
    '''
    return np.floor(round((val - min_val) / bin_size, 8)).astype('int')

def num_bins(min_val, max_val, bin_size):
    """
    Find number of bins in given range for a given bin size.
    """
    return np.floor(round((max_val - min_val) / bin_size, 8)).astype('int')

def scale(data):
    """
    Scale data to probability.
    """
    data = np.array(data, dtype = 'float')
    return (data - min(data)) / (np.sum(data) - min(data))

def get_hist(data, min_val, max_val, bin_size):
    """
    Generate histogram.
    """
    n_bins = num_bins(min_val, max_val, bin_size)
    hist = np.histogram(data, bins = n_bins,
                        range = (min_val, max_val))
    y = hist[0]
    x = hist[1]
    x = x + ((x[-1] - x[0]) / (len(x) - 1)) / 2
    x = x[np.where(x < x[-1])]
    return x, y
