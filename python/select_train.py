#! /Users/Rowen/Documents/Library/anaconda/bin/python

import sys
import numpy as np

# READ DATA

input_file = sys.argv[1]
data = np.genfromtxt(input_file, dtype = 'float', unpack = True)

# SAMPLE FRACTIONS

train_num = int(float(sys.argv[2]) * float(data.shape[1]))
valid_num = int(float(sys.argv[3]) * float(data.shape[1])) + train_num

# EXTRACT SAMPLES

index = np.arange(0, data.shape[1], 1)
np.random.shuffle(index)

train_data = data[:, index[:train_num]].transpose()
valid_data = data[:, index[train_num:valid_num]].transpose()
test_data = data[:, index[valid_num:]].transpose()

# OUTPUT DATA

train_file = input_file + '_train.dat'
valid_file = input_file + '_valid.dat'
test_file = input_file + '_test.dat'

np.savetxt(train_file, train_data, fmt = '%f')
np.savetxt(valid_file, valid_data, fmt = '%f')
np.savetxt(test_file, test_data, fmt = '%f')

print 'Saving training set to:', train_file
print 'Saving validation set to:', valid_file
print 'Saving testing set to:', test_file

