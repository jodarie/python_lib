#! /Users/Rowen/Documents/Library/anaconda/bin/python

# Assign Errors to Existing Magnitudes

import argparse, os.path
import numpy as np
import scipy.special as ss

from functions import astro

# Arguments

# READ ARGUMENTS

parser = argparse.ArgumentParser('MAG_ERRORS OPTIONS:')
parser.add_argument('-i', '--input_file', dest = 'input_file',
                    help = 'Input file name.')
parser.add_argument('-c', '--columns', dest = 'columns', type = int,
                  nargs = '+', help = 'File column numbers.')
parser.add_argument('-m', '--mag_limits', dest = 'mag_limits', type = float,
                  nargs = '+', help = 'Magnitude limits.')
parser.add_argument('-z', '--zero_pts', dest = 'zero_pts', type = float,
                  nargs = '+', help = 'Zero points.')
parser.add_argument('-s', '--survey', dest = 'survey', help = 'Survey name.')
opts = parser.parse_args()

if not opts.input_file:
    print 'Error: Input file not specified.'
    exit()
 
if not os.path.exists(opts.input_file):
    print 'Error: Input file not found.'
    exit()

if opts.columns:
    opts.columns = np.array(opts.columns, dtype = 'int') - 1
else:
    print 'Error: Column numbers not specified.'
    exit()

if opts.mag_limits:
    if len(opts.mag_limits) != len(opts.columns):
        print 'Error: Number of magnitude limits must',
        print 'match number of input columns.'
        exit()
        
if opts.zero_pts:
    if len(opts.zero_pts) != len(opts.columns):
        print 'Error: Number of zero points must',
        print 'match number of input columns.'
        exit()

# Survey Properties

def survey_assign(survey):
    if survey == 'sdss':
        mag_limits = np.array([22.0, 22.2, 22.2, 21.3, 20.5])
        zero_pts = np.array([24.63, 25.11, 24.80, 24.36, 22.83])
    if survey == 'des':
        mag_limits = np.array([25.2, 24.8, 24.0, 23.4, 21.7])
        zero_pts = np.array([25.296, 25.374, 25.379, 25.064, 24.154])
    if survey == 'euclid':
        mag_limits = np.array([24, 24, 24])
        zero_pts = np.array([24.25, 24.29, 24.92])
    return mag_limits, zero_pts

if opts.survey:
    opts.mag_limits, opts.zero_pts = survey_assign(opts.survey)

# Mag Errors

def capak_error(mag, mag_limit, zero_pt):
    flux_lim = astro.mag2flux(mag_limit, zero_pt)
    hold = 2.0 * np.random.random() - 1.0
    if hold > 0:
      hold = ss.erfinv(hold) * flux_lim
    else:
      hold = -1.0 * ss.erfinv(-1.0 * hold) * flux_lim
    flux = astro.mag2flux(mag, zero_pt) + hold
    fluxerror = (flux + flux_lim ** 2.0) ** 0.5 / 10.0
    obs_mag = astro.flux2mag(flux, zero_pt)
    mag_err = (2.5 / np.log(10.0)) * (fluxerror / flux)
    if flux < 0 or mag_err > 2.0:
      mag_err = mag_limit
    return obs_mag, mag_err

def add_errors(mags, mag_limits, zero_pts):
    new_mags = []
    new_mag_errs = []
    for i in range(len(mags)):
        new_mag, new_mag_err = capak_error(mags[i], mag_limits[i], zero_pts[i])
        new_mags.append(new_mag)
        new_mag_errs.append(new_mag_err)
    new_mags = np.array(new_mags)
    new_mag_errs = np.array(new_mag_errs)
    return new_mags, new_mag_errs

# Read Line

with open(opts.input_file) as f:
    for line in f:
        mags = np.array(line.split(' '), dtype = 'float')[opts.columns]
        new_mags, new_mag_errs = add_errors(mags, opts.mag_limits, opts.zero_pts)
        print ' '.join(map(str, mags)), ' '.join(map(str, new_mags)), \
          ' '.join(map(str, new_mag_errs))
        
