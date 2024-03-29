## @file nfw.py
#
#  NFW FUNCTIONS 
#
#  Functions for calculating NFW parameters.
#
#  REFERENCES:
#  1) Lokas and Mamon, Properties of spherical galaxies and clusters
#  with an NFW density profile, 2001, MNRAS 321, 155. (LM2001)
#  2) Mamon, Biviano and Boue, MAMPOSSt: Modelling Anisotropy and Mass
#  Profiles of Observed Spherical Systems - I. Gaussian 3D velocities,
#  2013, MNRAS, 429, 3079. (MBB2013)
#
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np

##
#  Function calculates the NFW profile
#  projected surface density.
#
#  Equations 42 and 43 from LM2001.
#
#  @param[in] t: t = R / r_s.
#
def nfw_proj_sd(t):

    if np.around(t, 8) == 0.0:
        sd = 0.0
        
    elif 1.0 - np.around(t, 8) == 0.0:
        sd = 1.0 / 3.0
        
    else:
        cm1 = 1.0
        if t < 1:
                cm1 = np.arccosh(1.0 / t)
        elif t > 1:
                cm1 = np.arccos(1.0 / t)
        sd = (1.0 - cm1 / np.sqrt(np.abs(t ** 2 - 1.0))) / (t ** 2 - 1.0)
    
    return sd / (2.0 * np.log(2.0) - 1.0)

##
#  Function calculates the NFW profile
#  projected mass along a cylinder in
#  units of M(r-2).
#
#  Equation B1 from MBB2013.
#
#  @param[in] t: t = R / r_s.
#
def nfw_proj_mass(t):

    if np.around(t, 8) == 0.0:
        mp = 0.0

    elif 1.0 - np.around(t, 8) == 0.0:
        mp = 1.0 - np.log(2.0)
        
    else:
        capc = 1.0
        if t > 1:
            capc = np.arccos(1.0 / t)
        elif t < 1:
            capc = np.arccosh(1.0 / t)
        mp = capc / np.sqrt(np.abs(t ** 2 - 1.0)) + np.log(t / 2.0)

    return mp / (np.log(2.0) - 0.5)

##
#  NFW number density for given projected
#  radius, scale radius and background
#  density.
#
#  @param[in] R_proj: Projected radius.
#  @param[in] r_s: Scale radius.
#
def nfw_num_density(R_proj, r_scale, bg_density):

    t_up = max(R_proj) / r_scale
    t_low = min(R_proj) / r_scale

    nfw_num = len(R_proj) - np.pi * r_scale ** 2 * (t_up ** 2 - t_low ** 2) * bg_density
    nfw_den = np.pi * r_scale ** 2 * (nfw_proj_mass(t_up) - nfw_proj_mass(t_low))

    if nfw_num <= 0:
        nfw_num = 0.1
    
    return nfw_num / nfw_den

##
#  Maximum liklihood for projected NFW
#  with background.
#
#  @param[in] r_scale_bg: Scale radius
#  and background density.
#  @param[in] R_proj: Projected radius.
#  @param[in] weights: Optional weights.
#
def nfw_proj_maxlik_bg(r_scale_bg, R_proj, weights = None):
    
    nfw_sd = np.array(map(nfw_proj_sd, R_proj / r_scale_bg[0]))

    nfw_nden = nfw_num_density(R_proj, r_scale_bg[0], r_scale_bg[1])
    
    prob = nfw_nden * nfw_sd + r_scale_bg[1]

    if weights is None:
        return np.sum(-np.log(prob))

    else:
        return np.sum(-np.log(prob) * weights)

##
#  Minimisation of projected NFW with background
#  density using adaptive grid search.
#
#  @param[in] R_proj: Projected radius.
#  @param[in] r_s: Scale radius.
#  @param[in] bg_density: Background density.
#  @param[in] weights: Optional weights.
#
def min_nfw_proj(R_proj, r_scale, bg_density = 0.0, fit_bg = True, weights = None):

    k = 0
    dist = 0.03
    r_scale_temp, bg_density_temp = r_scale, bg_density

    def get_bins(i, x):
        n_bins = 20 
        bin_step = 1.7 / n_bins / i ** 2 
        return 10.0 ** (np.log10(x) + (np.arange(n_bins) - (n_bins / 2.0)) * bin_step)

    def get_dist(best, temp):
        return np.sqrt(np.sum((best / temp - 1.0) ** 2))

    while k < 10 and dist > 0.02: 
        
        k += 1

        a_bins = get_bins(k, r_scale_temp)
        if fit_bg:
            b_bins = get_bins(k, bg_density_temp)
            
        ml_best = nfw_proj_maxlik_bg([r_scale_temp, bg_density_temp], R_proj, weights)
        r_scale_best, bg_density_best = r_scale_temp, bg_density_temp
        
        for i in range(len(a_bins)):
            if fit_bg:
                for j in range(len(b_bins)):
                    mli = nfw_proj_maxlik_bg([a_bins[i], b_bins[j]], R_proj, weights)
                    if mli < ml_best: 
                        ml_best = mli
                        r_scale_best, bg_density_best = a_bins[i], b_bins[j]
            else:
                mli = nfw_proj_maxlik_bg([a_bins[i], bg_density_temp], R_proj, weights)
                if mli < ml_best:
                    ml_best = mli
                    r_scale_best = a_bins[i]

        if fit_bg:
            dist = get_dist(np.array([r_scale_best, bg_density_best]),
                            np.array([r_scale_temp, bg_density_temp]))
        else:
            dist = get_dist(r_scale_best, r_scale_temp)
            
        print '', k, ': Output Solution', [r_scale_best, bg_density_best]

        r_scale_temp, bg_density_temp = r_scale_best, bg_density_best

    print 'ML-Best =', ml_best
    
    return r_scale_best, bg_density_best
