#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np

data = np.genfromtxt(sys.argv[1], dtype="float", unpack = True)
z = data[0,:]
mass = data[1,:]
rich = np.log10(data[2,:])

def y(x, m, c):
    return m * x + c

for i in range(len(z)):
    print z[i], mass[i], 10 ** (rich[i] / y(z[i] ,float(sys.argv[2]), float(sys.argv[3])))

