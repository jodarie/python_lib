#! /Users/Rowen/Documents/Library/anaconda/bin/python

import numpy as np
import sys

data = np.genfromtxt(sys.argv[1], unpack = True, dtype = 'float')

column = int(sys.argv[2]) - 1

print 'Min = ', np.min(data[column, :])
print 'Max = ', np.max(data[column, :])
