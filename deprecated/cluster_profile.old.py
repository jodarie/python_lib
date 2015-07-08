#! /Users/Rowen/Documents/Library/anaconda/bin/python

import sys
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt

from functions import astro, cosmo2
from biviano_py import profmaxliknfwp

# adopted cosmology
h0 = 70.0 
omega_m = 0.3
omega_l = 0.7 

# read cluster file
data = np.genfromtxt(sys.argv[1], unpack = True,
                             dtype = 'float')
cluster = np.genfromtxt(sys.argv[2], unpack = True,
                             dtype = 'float')
cluster_num = float(sys.argv[3])
sncl = data[8, data[0] == cluster_num]

index = cluster[0] == cluster_num

# cluster properties
racl = np.median(cluster[4, index])
decl = np.median(cluster[5, index])
zcl = np.median(cluster[6, index])
dists = astro.ang_sep([cluster[4, index], cluster[5, index]], [racl,  decl]) * 60.0
area = np.pi * np.mean(dists) ** 2 # arcmin^2

# Mpc per arcmin scale
da = cosmo2.d_angdi(zcl, omega_m, omega_l) * cosmo2.d_H(h0)
mpca = astro.deg2rad(da / 60.0)

# coordinates wrt cluster center
gx = 60. * (racl - cluster[4, index]) * np.cos(astro.deg2rad(cluster[5, index]))
gy = 60. * (cluster[5, index] - decl)   # arcmin      

# adaptive kernel density map
gadk = ss.gaussian_kde(np.vstack([gx * mpca, gy * mpca]), bw_method = 'silverman')

xedges = np.linspace(-0.6, 0.6, 100)
yedges = np.linspace(-1.0, 1.0, 100)
xx, yy = np.meshgrid(xedges, yedges)
gridpoints = np.array([xx.ravel(), yy.ravel()])
zz = np.reshape(gadk(gridpoints), xx.shape)
indices = np.unravel_index(zz.argmax(), zz.shape)

# distance in Mpc from cluster center
adkcx = xedges[indices[1]]
adkcy = yedges[indices[0]]

#adkcx = 0.0121250 #HOT WIRED#   
#adkcy = 0.0900006 #HOT WIRED#   

rmpc = mpca * np.sqrt((gx - adkcx) ** 2 + (gy - adkcy) ** 2)
cm = np.ones(len(cluster[0]))

# plot
plt.figure(1)
plt.plot(gx * mpca, gy * mpca, 'k+')
plt.xlabel('x [Mpc]')
plt.ylabel('y [Mpc]')
plt.xlim(-0.6, 0.6)
plt.ylim(-1.0, 1.0)
plt.contourf(xedges, yedges, zz, 20)
plt.plot(adkcx, adkcy, 'wx')
plt.colorbar()

# need to estimate a guess for the background density;
# we use the S/N, the number of galaxies assigned to
# the cluster within the cluster area. We write:
# S/N = Ncl / sqrt(Nfd), assuming Poissonian error for
# fd counts. Then we assume that Ncl = number of galaxies
# within the area.
rlim = np.sqrt(area / np.pi) * mpca

w = np.where(rmpc <= rlim)[0]
narea = len(w)

nfd = (narea / sncl) ** 2

# free parameters in the fit: scale radius and bkd number density (in
# Mpc and Mpc^-2 units)
parinp = [0.3, nfd / (area * mpca ** 2)] 
npfree = 2

rnu0, bkd0, chiq, chiqp, rd, dp, edp, xf, yf = profmaxliknfwp.profmaxliknfwp(rmpc[w], cm[w], 1, parinp, npfree, 100)

print rd, dp, edp, rnu0, bkd0

# plot
plt.figure(2)
plt.plot(rd, dp, 'bx')
plt.errorbar(rd, dp, yerr = edp, linestyle = 'None')
plt.xlim(np.min(rd) / 3., 1.5 * np.max(rd))
plt.ylim(np.min(dp) / 2., np.max(dp) * 2.)
plt.text(0.3 * max(rd), max(dp), r'$r_s =$ ' + str(rnu0))
plt.text(0.3 * max(rd), 0.77 * max(dp), r'$Bkd =$ ' + str(bkd0))
plt.xlabel('R [Mpc]')
plt.ylabel('number density [Mpc$^{-2}$]')
plt.xscale('log')
plt.yscale('log')
plt.title('projected NFW best-fit')
for i in range(len(rd)):
    plt.plot(rd[i], rd[i]),
plt.plot(xf, yf, 'r--')
plt.show()
