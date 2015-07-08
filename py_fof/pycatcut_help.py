#################
# PYCATCUT HELP #
#################

def print_help():
    '''
    Function that prints the help information for pycatcut.py.
    '''
    print '===================='
    print ' PYCATCUT V1.0 2014 '
    print '===================='
    print ''
    print '1) To run the code:'
    print ''
    print '  >> python pycatcut -f <FILE_NAME> --ra_col <RA_COLUMN> --dec_col <DEC_COLUMN> --ra_bin <NUM_RA_BINS> --dec_bin <NUM_DEC_BINS>'
    print ''
    print '2) Options:'
    print ''
    print '  --output      (set output file prefix)'
    print '  --ra_col      (set column number of RA in the catalogue)'
    print '  --dec_col     (set column number of DEC in the catalogue)'
    print '  --ra_lower    (set lower limit in RA range) [min RA]'
    print '  --ra_upper    (set upper limit in RA range) [max RA]'
    print '  --dec_lower   (set lower limit in DEC range) [min DEC]'
    print '  --dec_upper   (set upper limit in DEC range) [max DEC]'
    print '  --ra_bin      (set number of RA bins)'
    print '  --dec_bin     (set number of DEC bins)'
    print '  --ra_overlap  (set overlap between RA bins) [default 1.0]'
    print '  --dec_overlap (set overlap between DEC bins) [default 1.0]'
    print ''
    print '3) Examples:'
    print ''
    print '    >> python pycatcut -f 2slaq_spec.dat --ra_col 2 --dec_col 3 --ra_bin 3 --dec_bin 3'
    print ''
    exit()

def check_errors(ra_col, dec_col, ra_bin, dec_bin):
    '''
    Function that checks for errors in the arguments.
    '''
    if ra_col == -1:
        print 'ERROR! RA column number not specified.'
        exit()
    elif dec_col == -1:
        print 'ERROR! DEC column number not specified.'
        exit()
    elif ra_col == dec_col:
        print 'ERROR! Same column number give for RA and DEC fields.'
        exit()
    elif ra_bin == -1 or ra_bin < 1:
        print 'ERROR! Number of RA bins not specified or invalid.'
        exit()
    elif dec_bin == -1 or dec_bin < 1:
        print 'ERROR! Number of DEC bins not specified or invalid.'
        exit()
