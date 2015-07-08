#! /Users/Rowen/Documents/Library/anaconda/bin/python

################
# ANALYZ V.1.0 #
################

import math
import numpy as np
import numpy.random as npr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import sys, getopt

def usage():
   print 'analyz.py -i <input_file> -o <output_file> -x <x_column> -y <y_column> -b <bin_size>'

def percentile(N, percent): 
    N=np.sort(N)
    k = (len(N)-1) * percent
    f = np.floor(k)
    c = np.ceil(k)
    if f == c: return N[int(k)]
    d0 = N[int(f)] * (c-k)
    d1 = N[int(c)] * (k-f)
    return d0+d1

def std68(N): 
    calc68 = 0.5*(percentile(N,0.84)-percentile(N,0.16))
    return calc68
 
def outfrac(data,limit):
   value = ((data>limit).sum(0)+(data<-limit).sum(0))/float(len(data))
   return value

def bootstrap(data,num_samples,statistic,alpha,limit):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    n = len(data)
    if statistic==outfrac:
       thestat = statistic(data,limit)
    else:
       thestat = statistic(data)
    idx = npr.randint(0,n,(num_samples,n))
    samples = data[idx]
    if statistic==std68:
       samples = np.transpose(samples)
       stat = np.sort(statistic(samples))
    elif statistic==outfrac:
       samples = np.transpose(samples)
       stat = np.sort(statistic(samples,limit))
    else: stat = np.sort(statistic(samples, 1))
    low = stat[int((alpha/2.0)*num_samples)]
    high = stat[int((1-alpha/2.0)*num_samples)]
    value1 = thestat-low
    value2 = high-thestat
    return np.mean([value1,value2]) 
   
def main(argv):
   
   ##################
   # READ ARGUMENTS #
   ##################
   
   try:
      opts, args = getopt.getopt(argv,"i:o:x:y:b:h",["ifile=","ofile=","zspec=","zphot=","binsize=","help"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   if len(opts) < 4:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         usage()
         sys.exit(2)
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-x", "--zspec"):
         xcol = arg
      elif opt in ("-y", "--zphot"):
         ycol = arg
      elif opt in ("-b", "--binsize"):
         binsize = arg
      else:
         usage()
         sys.exit(2)

   xcol = int(xcol) - 1
   ycol = int(ycol) - 1
   binsize = float(binsize)
   
   #############
   # READ DATA #
   #############
         
   data = np.genfromtxt(inputfile, unpack = True)

   z_spec = data[xcol,:]
   z_phot = data[ycol,:]

   upperlimit = max(z_phot)
   lowerlimit = min(z_phot)
   upperlimit2 = max(z_spec)+0.1
   lowerlimit2 = min(z_spec)
   xyline = (lowerlimit,upperlimit)
   xyline2 = (lowerlimit2,upperlimit2)

   ##############
   # STATISTICS #
   ##############

   dz1 = (z_spec-z_phot)
   dz2 = (z_spec-z_phot)/(1+z_spec)

   #bootsrap properties
   num_samples = 1000
   confidence = 0.05

   #bias
   bias = np.mean(dz1)
   biaserr = bootstrap(dz1,num_samples,np.mean,confidence,0)   
   dbias = np.mean(dz2)
   dbiaserr = bootstrap(dz2,num_samples,np.mean,confidence,0)

   #sigma
   sigma = np.std(dz1)
   sigmaerr = bootstrap(dz1,num_samples,np.std,confidence,0)
   dsigma = np.std(dz2)
   dsigmaerr = bootstrap(dz2,num_samples,np.std,confidence,0)

   #sigma68
   sigma68 = std68(dz1)   
   sigma68err = bootstrap(dz1,num_samples,std68,confidence,0)
   dsigma68 = std68(dz2)
   dsigma68err = bootstrap(dz2,num_samples,std68,confidence,0)

   #outliers
   outlimit = 2*sigma
   doutlimit = 2*dsigma
   outfraction = outfrac(dz1,outlimit)
   outfractionerr = bootstrap(dz1,num_samples,outfrac,confidence,outlimit)
   
   print "========================================="
   print "Average Bias: %(bval)0.5f +/- %(beval)0.5f" % {'bval': bias, 'beval': biaserr}
   print "Average Bias/(1+z): %(bval)0.5f +/- %(beval)0.5f" % {'bval': dbias, 'beval': dbiaserr}
   print "Sigma: %(bval)0.5f +/- %(beval)0.5f" % {'bval': sigma, 'beval': sigmaerr}
   print "Sigma/(1+z): %(bval)0.5f +/- %(beval)0.5f" % {'bval': dsigma, 'beval': dsigmaerr}
   print "Sigma68: %(bval)0.5f +/- %(beval)0.5f" % {'bval': sigma68, 'beval': sigma68err}
   print "Sigma68/(1+z): %(bval)0.5f +/- %(beval)0.5f" % {'bval': dsigma68, 'beval': dsigma68err}
   print "Outlier Fraction: %(bval)0.5f +/- %(beval)0.5f" % {'bval': outfraction, 'beval': outfractionerr}
   print "========================================="

   #redshift bins
   nbins = int(math.floor((max(z_spec)-min(z_spec))/binsize))+1

   #statistics as a function of z
   xvals = [] # array with z_spec bin centres
   binmin = [] # array with z_spec bin lower limit
   binmax = [] # array with z_spec bin upper limit
   data2d = [] # 2D array with binned data
   ddata2d = [] # 2D array with binned data
   bias2d = [] # array bias as a function of z
   bias2derr = [] # array bias error as a function of z
   dbias2d = [] # array bias/(1+z) as a function of z
   dbias2derr = [] # array bias/(1+z) error as a function of z
   sigma2d = [] # array sigma as a function of z
   sigma2derr = [] # array sigma error as a function of z
   dsigma2d = [] # array sigma/(1+z) as a function of z
   dsigma2derr = [] # array sigma/(1+z) error as a function of z
   sigma682d = []  # array sigma_68 as a function of z
   sigma682derr = [] # array sigma_68 error as a function of z
   dsigma682d = [] # array sigma_68/(1+z) as a function of z
   dsigma682derr = [] # array sigma_68/(1+z) error as a function of z
   outfraction2d = [] # array outlier fraction as a function of z
   outfraction2derr = [] # array outlier fraction error as a function of z
   
   for i in range(nbins):
      xvals.append(min(z_spec)+binsize*(i+0.5))
      binmin.append(min(z_spec)+i*binsize)
      binmax.append(min(z_spec)+(i+1)*binsize)
      data2d.append([])
      ddata2d.append([])

   xvals = np.array(xvals)
   binmin = np.array(binmin)
   binmax = np.array(binmax)
      
   for i in range(len(dz1)): 
      for j in range(nbins):
         if z_spec[i]>=binmin[j] and z_spec[i]<binmax[j]:
            data2d[j].append(dz1[i])
            ddata2d[j].append(dz2[i])

   data2d = np.array(data2d)
   ddata2d = np.array(ddata2d)
   for i in range(nbins):
      data2d[i] = np.array(data2d[i])
      ddata2d[i] = np.array(ddata2d[i])

   for i in range(nbins):
      bias2d.append(np.mean(data2d[i]))
      bias2derr.append(bootstrap(data2d[i],num_samples,np.mean,confidence,0))
      dbias2d.append(np.mean(ddata2d[i]))
      dbias2derr.append(bootstrap(ddata2d[i],num_samples,np.mean,confidence,0))
      sigma2d.append(np.std(data2d[i]))
      sigma2derr.append(bootstrap(data2d[i],num_samples,np.std,confidence,0))
      dsigma2d.append(np.std(ddata2d[i]))
      dsigma2derr.append(bootstrap(ddata2d[i],num_samples,np.std,confidence,0))
      sigma682d.append(std68(data2d[i]))
      sigma682derr.append(bootstrap(data2d[i],num_samples,std68,confidence,0))
      dsigma682d.append(std68(ddata2d[i]))
      dsigma682derr.append(bootstrap(ddata2d[i],num_samples,std68,confidence,0))
      outfraction2d.append(outfrac(data2d[i],outlimit))
      outfraction2derr.append(bootstrap(data2d[i],num_samples,outfrac,confidence,outlimit))
      
   bias2d = np.array(bias2d)   
   bias2derr = np.array(bias2derr)
   dbias2d = np.array(dbias2d)   
   dbias2derr = np.array(dbias2derr)
   sigma = np.array(sigma2d)   
   sigmaerr = np.array(sigma2derr)
   dsigma = np.array(dsigma2d)   
   dsigmaerr = np.array(dsigma2derr)
   sigma68 = np.array(sigma682d)   
   sigma68err = np.array(sigma682derr)
   dsigma68 = np.array(dsigma682d)   
   dsigma68err = np.array(dsigma682derr)
   outfraction2d = np.array(outfraction2d)
   outfraction2derr = np.array(outfraction2derr)

   print outfraction2derr
   
   ################
   # OUTPUT PLOTS #
   ################

   print "Making plots..."

   #Scatter Plot

   fig1 = plt.figure() 
   fig1.subplots_adjust(hspace = .001)
   gs1 = gridspec.GridSpec(3, 1, height_ratios=[3,1,1])
   ax1 = fig1.add_subplot(gs1[0])
   ax2 = fig1.add_subplot(gs1[1],sharex=ax1)
   ax3 = fig1.add_subplot(gs1[2],sharex=ax1)

   ax1.scatter(z_spec, z_phot, s=2, c='0', lw=0)
   ax1.plot(xyline,xyline,'--',c='r',linewidth=1.5)
   ax1.set_xlim(lowerlimit,upperlimit)
   ax1.set_ylim(lowerlimit,upperlimit)
   ax1.yaxis.set_major_locator(MaxNLocator(6))
   ax1.set_title('Photo-z Scatter')
   ax1.set_ylabel('$z_{phot}$')
   ax1.xaxis.set_visible(False)

   ax2.scatter(z_spec, dz1, s=2, c='0', lw=0)
   ax2.plot(xyline,(0,0),'-.',c='b',linewidth=1.5)
   ax2.plot(xyline,(bias,bias),'--',c='r',linewidth=1.5)
   ax2.axhspan(-outlimit,outlimit,facecolor='b',alpha=0.2,edgecolor='none')
   ax2.set_xlim(lowerlimit,upperlimit)
   ax2.set_ylim(min(dz1),max(dz1))
   ax2.yaxis.set_major_locator(MaxNLocator(5))
   ax2.set_ylabel('$\Delta z$')
   ax2.xaxis.set_visible(False)

   ax3.scatter(z_spec, dz2, s=2, c='0', lw=0)
   ax3.plot(xyline,(0,0),'-.',c='b',linewidth=1.5)
   ax3.plot(xyline,(bias,bias),'--',c='r',linewidth=1.5)
   ax3.axhspan(-doutlimit,doutlimit,facecolor='b',alpha=0.2,edgecolor='none')
   ax3.set_xlim(lowerlimit,upperlimit)
   ax3.set_ylim(min(dz2),max(dz2))
   ax3.xaxis.set_major_locator(MaxNLocator(6))
   ax3.yaxis.set_major_locator(MaxNLocator(5))
   ax3.set_xlabel('$z_{spec}$')
   ax3.set_ylabel('$\Delta z/(1+z)$')

   fig1_name = outputfile + "_scatter"
   fig1.savefig(fig1_name)
   print "Output saved to", fig1_name

   #Density Plot
   
   heat, xedges, yedges = np.histogram2d(z_spec, z_phot, bins=(50))
   extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

   fig2 = plt.figure() 
   ax4 = fig2.add_axes([0.1,0.1,0.6,0.8])
   im = ax4.imshow(heat, extent=extent, interpolation='spline16', cmap=cm.jet)
   ax4.set_xlim(lowerlimit,upperlimit)
   ax4.set_ylim(lowerlimit,upperlimit)
   ax4.set_title('Photo-z Scatter Density')
   ax4.set_xlabel('$z_{spec}$')
   ax4.set_ylabel('$z_{phot}$')
   cax = fig2.add_axes([0.7, 0.1, 0.05, 0.8])
   cbar = fig2.colorbar(im,cax=cax)
   cbar.set_label('Photoz Density')

   fig2_name = outputfile + "_density"
   fig2.savefig(fig2_name)
   print "Output saved to", fig2_name

   #N(z) Plot

   nz_binsize = 0.01
   nz_nbins1 = (max(z_spec)-min(z_spec))/nz_binsize
   nz_nbins2 = (max(z_phot)-min(z_phot))/nz_binsize

   fig3 = plt.figure() 
   plt.hist(z_spec, nz_nbins1, normed=1, facecolor='green', alpha=0.5, label='$z_{spec}$')
   plt.hist(z_phot, nz_nbins2, normed=1, facecolor='red', alpha=0.5, label='$z_{phot}$')
   plt.legend()
   plt.xlabel('$z$')
   plt.title('N(z) Distribution')

   fig3_name = outputfile + "_nz"
   fig3.savefig(fig3_name)
   print "Output saved to", fig3_name

   #bias(z)/sigama(z)/sigma68(z)/Outliers(z) Plots

   fig4 = plt.figure() 
   fig4.subplots_adjust(hspace = .001)
   plt.rc('legend',**{'fontsize':6})
   
   gs2 = gridspec.GridSpec(4, 1, height_ratios=[1,1,1,1])
   ax5 = fig4.add_subplot(gs2[0])
   ax6 = fig4.add_subplot(gs2[1],sharex=ax5)
   ax7 = fig4.add_subplot(gs2[2],sharex=ax5)
   ax8 = fig4.add_subplot(gs2[3],sharex=ax5)
   
   ax5.plot(xvals,bias2d,linewidth=1.5,c='b',label='Bias')
   ax5.errorbar(xvals,bias2d,yerr=bias2derr,fmt='bx')
   ax5.plot(xvals,dbias2d,linewidth=1.5,c='g',label='Bias/$(1+z)$')
   ax5.errorbar(xvals,dbias2d,yerr=dbias2derr,fmt='gx')
   ax5.plot(xyline2,(0,0),'--',c='r',linewidth=1.5)
   ax5.set_xlim(lowerlimit2,upperlimit2)
   ax5.set_ylabel('Bias')
   ax5.set_title('Metrics(z)')
   ax5.legend()
   ax5.yaxis.set_major_locator(MaxNLocator(4))
   ax5.xaxis.set_visible(False)

   ax6.plot(xvals,sigma2d,linewidth=1.5,c='b',label='$\sigma$')
   ax6.errorbar(xvals,sigma2d,yerr=sigma2derr,fmt='bx')
   ax6.plot(xvals,dsigma2d,linewidth=1.5,c='g',label='$\sigma$/(1+z)')
   ax6.errorbar(xvals,dsigma2d,yerr=dsigma2derr,fmt='gx')
   ax6.plot(xyline2,(0,0),'--',c='r',linewidth=1.5)
   ax6.set_xlim(lowerlimit2,upperlimit2)
   ax6.set_ylabel('$\sigma$')
   ax6.legend()
   ax6.yaxis.set_major_locator(MaxNLocator(3))
   ax6.xaxis.set_visible(False)

   ax7.plot(xvals,sigma682d,linewidth=1.5,c='b',label='$\sigma_{68}$')
   ax7.errorbar(xvals,sigma682d,yerr=sigma682derr,fmt='bx')
   ax7.plot(xvals,dsigma682d,linewidth=1.5,c='g',label='$\sigma_{86}/(1+z)$')
   ax7.errorbar(xvals,dsigma682d,yerr=dsigma682derr,fmt='gx')
   ax7.plot(xyline2,(0,0),'--',c='r',linewidth=1.5)
   ax7.set_xlim(lowerlimit2,upperlimit2)
   ax7.set_ylabel('$\sigma_{68}$')
   ax7.legend()
   ax7.yaxis.set_major_locator(MaxNLocator(3))
   ax7.xaxis.set_visible(False)

   ax8.plot(xvals,outfraction2d,linewidth=1.5,c='b')
   ax8.errorbar(xvals,outfraction2d,yerr=outfraction2derr,fmt='bx')
   ax8.plot(xyline2,(0,0),'--',c='r',linewidth=1.5)
   ax8.set_xlim(lowerlimit2,upperlimit2)
   ax8.set_ylabel('Outliers')
   ax8.set_xlabel('$z_{spec}$')
   ax8.yaxis.set_major_locator(MaxNLocator(3))

   fig4_name = outputfile + "_bias"
   fig4.savefig(fig4_name)
   print "Output saved to", fig4_name   

if __name__ == "__main__":
   main(sys.argv[1:])
