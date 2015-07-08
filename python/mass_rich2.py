#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np, matplotlib as mp, scipy as sp
from scipy import stats
from scipy.odr import *
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

# FIT MODEL

def f(a, b, x):
    return b * x + a

# READ FILE

file_name = sys.argv[1]

data = np.genfromtxt(file_name, dtype="float", unpack = True)

z = data[0,:]
mass = data[1,:]
rich_log = np.log10(data[2,:])

# LINEAR REGRESSION

b1, a1, r_value, p_value, std_err = stats.linregress(mass, rich_log)
b2, a2, r_value, p_value, std_err = stats.linregress(rich_log, mass)

# PLOT

x1 = [min(mass), max(mass)]
y1 = [f(a1, b1, x1[0]), f(a1, b1, x1[1])]

y2 = [min(rich_log), max(rich_log)]
x2 = [f(a2, b2, y2[0]), f(a2, b2, y2[1])]

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.scatter(mass, rich_log, s=2, c='0', lw=0)
beta = (r'$\beta$')
ax.plot(x1, y1, label = 'P(Y|X), ' + beta + ' = ' + str(round(b1, 3)))
ax.plot(x2, y2, label = 'P(X|Y), ' + beta + ' = ' + str(round(1.0 / b2, 3)))
ax.set_xlabel('Log$_{10}$ M$_{Halo}$')
ax.set_ylabel('Log$_{10}$ N$_{gal}$')
ax.set_title(file_name + " mass-richness relation")
ax.legend(loc = 2, prop = {'size':10})
fig_name = file_name + '_mass_v_rich_plot.1.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name
