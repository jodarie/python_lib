#! /Users/Rowen/Documents/Library/anaconda/bin/python

import sys
import numpy as np
import pandas as pd

file_name = sys.argv[1]

test_data = np.genfromtxt(file_name, unpack = True, dtype = 'S')[[3, 4, 5, 6]]

test_data2 = pd.DataFrame({'id' : test_data[0], 
                           'ra' : np.array(test_data[1], dtype = 'float64'),
                           'dec' : np.array(test_data[2], dtype = 'float64'),
                           'z' : np.array(test_data[3], dtype = 'float64')})


print test_data2.mean()[2]
