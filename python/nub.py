#! /Users/Rowen/Documents/Library/anaconda/bin/python

import numpy as np
import sys

file_name = sys.argv[1]

data = np.genfromtxt(file_name, unpack = True, dtype="S")

count = np.array(data[0,:], dtype = 'int')
ids = data[1,:]

for i in range(len(ids)):
    for j in range(count[i]):
        print ids[i], count[i]
