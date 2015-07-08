#! /Users/Rowen/Documents/Library/anaconda/bin/python

import math, sys
import numpy as np, matplotlib as mp, scipy as sp
from scipy import stats
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from functions import fits
from biviano_py import biweight
from web_functions import robust

# ERRORS

def get_err(x, y):
    xerr = 1.0 / (np.sqrt(x) * np.log(10))
    yerr = np.zeros(len(y)) + 0.01
    return xerr, yerr

# READ FILE

file_name = sys.argv[1]

data = np.genfromtxt(file_name, dtype="float", unpack = True)

z = data[0,:]
mass = data[1,:]
rich = data[2,:]
rich_log = (np.log10(rich) - min(np.log10(rich))) / (max(np.log10(rich)) - min(np.log10(rich)))

rich_log_err, mass_err = get_err(rich, mass)

# MASS LIMIT

mass_limit = float(sys.argv[2])
mass_index = (mass >= mass_limit)

# PLOT COLOUR

plot_colour = sys.argv[3]

# ODR FIT

beta = fits.fit_odr(rich_log[mass_index], mass[mass_index], rich_log_err[mass_index], mass_err[mass_index], fits.linear_fit)

# SCATTER

m_step = 0.25
m_bins = min(mass) + m_step * np.arange(10, dtype = 'float')
r_bins = (m_bins - beta[1]) / beta[0]

dev = mass - beta[1] - beta[0] * rich_log

scat1 = []
scat2 = []
scat3 = []

for i in range(len(r_bins) - 1):
    index = np.where(((rich_log >= r_bins[i]) & (rich_log < r_bins[i + 1])))[0]
    scat1.append(robust.std(dev[index]))
    scat2.append(np.std(dev[index]))
    scat3.append(biweight.bwtav(dev[index])[3])

scat1 = np.array(scat1)
scat2 = np.array(scat2)
scat3 = np.array(scat3)

m_bins = m_bins[:len(m_bins) - 1] + m_step / 2.0
    
# PLOT

x = [min(rich_log), max(rich_log)]
y = [fits.linear_fit(beta, x[0]), fits.linear_fit(beta, x[1])]
m_lim = [mass_limit, mass_limit]

## Rich v Mass

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.scatter(rich_log, mass, s=2, c='0', lw=0, color = plot_colour)
ax.set_ylim(min(mass), max(mass))
ax.set_xlabel('Log$_{10}$ Cluster Mass Proxy')
ax.set_ylabel('Log$_{10}$ Halo Mass')
fig_name = file_name + '_mass_v_rich_plot0_' + str(mass_limit) + '.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.scatter(rich_log, mass, s=2, c='0', lw=0)
ax.plot(x, y)
ax.plot(x, m_lim, '--r')
ax.set_ylim(min(mass), max(mass))
ax.set_xlabel('Log$_{10}$ Cluster Mass Proxy')
ax.set_ylabel('Log$_{10}$ Halo Mass')
ax.set_title(file_name + " richness-mass relation")
fig_name = file_name + '_mass_v_rich_plot1_' + str(mass_limit) + '.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

print "Fit: ", beta

## Mass v Scatter

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.plot(m_bins, scat1, label = 'Robust')
ax.plot(m_bins, scat2, label = 'Numpy')
ax.plot(m_bins, scat3, label = 'Biviano')
ax.set_xlim(13.5, 15.0)
ax.set_ylim(0.0, 0.4)
ax.legend()
fig_name = file_name + '_mass_v_scatter_plot_' + str(mass_limit) + '.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

################################################

mass_low_limit = 14.2
mass_high_limit = 14.4

z_mid = z[((mass>=mass_low_limit) & (mass<mass_high_limit))]
rich_mid = rich_log[((mass>=mass_low_limit) & (mass<mass_high_limit))]

z_bin_size = 0.5

n_z_bins = int(math.floor((max(z) - min(z)) / z_bin_size))

y_z = []
labels = []
z_list = []
rich_list = []
rich_zero = 0

for i in range(n_z_bins - 1):
    low = i * z_bin_size
    high = (i + 1) * z_bin_size
    labels.append(str(low) + ' <= z < ' + str(high))
    x_val = rich_log[((z >= low) & (z < high) & mass_index)]
    y_val = mass[((z >= low) & (z < high) & mass_index)]
    x_val_err = rich_log_err[((z >= low) & (z < high) & mass_index)]
    y_val_err = mass_err[((z >= low) & (z < high) & mass_index)]
    beta_val = fits.fit_odr(x_val, y_val, x_val_err, y_val_err, fits.linear_fit)
    y_z.append([fits.linear_fit(beta_val, x[0]), fits.linear_fit(beta_val, x[1])])
    rich_val = np.mean(rich_mid[((z_mid >= low) & (z_mid < high))])
    if i == 0:
        rich_zero = rich_val
    rich_list.append(rich_val / rich_zero)
    z_list.append(low + (z_bin_size / 2.))

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.scatter(rich_log, mass, s=2, c='0', lw=0)
for i in range(n_z_bins - 1):
   ax.plot(x, y_z[i], label = labels[i])
ax.set_ylim(min(mass), max(mass))
ax.set_xlabel('Log$_{10}$ Cluster Mass Proxy')
ax.set_ylabel('Log$_{10}$ Halo Mass')
ax.legend(loc=2,prop={'size':8})
fig_name = file_name + '_mass_v_rich_plot2_' + str(mass_limit) + '.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

coefficients = np.polyfit(z_list, rich_list, 1)
polynomial = np.poly1d(coefficients)
fit_vals = polynomial(z_list)

fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.plot(z_list, rich_list, label = "Average Ngal")
ax.plot(z_list, fit_vals, label = "Fit", color = 'red')
ax.set_xlabel('z')
ax.set_ylabel('Ngal')
ax.legend(loc=1,prop={'size':8})
fig_name = 'z_v_rich_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

print coefficients
