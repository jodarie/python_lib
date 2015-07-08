#! /Users/Rowen/Documents/Library/anaconda/bin/python

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser('QUICK-HIST OPTIONS:')
parser.add_argument('-i', '--input_files', dest = 'input_files',
                    nargs = '+', help = 'Input file name.')
parser.add_argument('-c', '--columns', dest = 'columns', type = int,
                    nargs = '+', help = 'File column number.')
parser.add_argument('-r', '--range', dest = 'range', default = [0.0, 1.0],
                    type = float, nargs = '+',
                    help = 'X-range. [Default: 0.0 1.0]')
parser.add_argument('-b', '--bins', dest = 'bins', default = 50,
                    type = int, help = 'Number of bins. [Default: 50]')
parser.add_argument('-t', '--type', dest = 'type', default = 'step',
                    help = 'Histogram type. [Default: \'bar\']')
parser.add_argument('-x', '--xlabel', dest = 'xlabel', default = 'Data',
                    help = 'X-axis label.')
opts = parser.parse_args()

if len(opts.input_files) == 1:
    data = np.genfromtxt(opts.input_files[0], unpack = True, dtype = 'float')
    for col in opts.columns:
        plt.hist(data[col - 1], bins = int(opts.bins), range = opts.range,
                histtype = opts.type, label = str(col))

else:
    for i in range(len(opts.input_files)):
        data = np.genfromtxt(opts.input_files[i], unpack = True, dtype = 'float')
        plt.hist(data[opts.columns[i] - 1], bins = int(opts.bins),
                 range = opts.range, histtype = opts.type, label = opts.input_files[i])
    
plt.xlim(opts.range)
plt.xlabel(opts.xlabel)
plt.legend()
plt.show()
