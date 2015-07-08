def confidence(r_scale, bg_density, R_proj, grid, dof):
    
    def get_range(value, grid):
        return tuple(10.0 ** (np.log10(value) + np.array([-grid / 2, grid / 2, -grid / 2]) * 1.4 / grid))

    ranges = (get_range(r_scale, grid), get_range(bg_density, grid))

    res = brute(nfw_proj_maxlik_bg, ranges, args = (R_proj, ), full_output = True)

    x, y = np.array(res[2]), np.array(res[3])
    x = x.reshape(x.shape[0], x.shape[1] * x.shape[2]).T
    y = y.flatten()

    index = (2.0 * (y - np.min(y)) <= chi2.isf(0.68, dof))

    lims = x[index, 0]

    def get_range2(value, grid):
        return 10.0 ** (np.log10(value) + np.arange(-grid / 2, grid / 2 + 1, 1) * 1.4 / grid)
    
    r_scale_bg_range = np.array([get_range2(r_scale, grid), get_range2(bg_density, grid)])
    r_scale_bg_range = np.array(list(product(*r_scale_bg_range)))
    
    gr = np.array(map(partial(nfw_proj_maxlik_bg, R_proj = R_proj), r_scale_bg_range))

    limits = r_scale_bg_range[(2.0 * (gr - np.min(gr)) <= chi2.isf(0.68, dof)), 0]

    print np.min(limits), np.max(limits), np.min(lims), np.max(lims)

    exit()
    
    return np.min(limits), np.max(limits)
