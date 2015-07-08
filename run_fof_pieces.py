#! /Users/Rowen/Documents/Library/anaconda/bin/python

import argparse
import subprocess as sp

# OPTIONS

parser = argparse.ArgumentParser('CODE OPTIONS:')
parser.add_argument('-r', '--link_r', dest = 'link_r',
                    type = float,
                    help = 'Angular linking paramter.')
parser.add_argument('-z', '--link_z', dest = 'link_z',
                    type = float,
                    help = 'L.O.S. linking paramter.')
parser.add_argument('-n', '--n_pieces', dest = 'n_pieces',
                    type = int,
                    help = 'Number of pieces.')
parser.add_argument('-d', '--dir', dest = 'dir', default = './',
                    help = 'Directory where pieces are stored.')
parser.add_argument('-m', '--mode', dest = 'fof_mode',
                    default = 'phot',
                    help = 'FoF Mode [Default: phot].')
parser.add_argument('-l', '--limits', dest = 'limits', nargs = '+',
                    default = [0.0, 3.0],
                    help = 'Redshift limits (min_z max_z).')
parser.add_argument('--min_ngal', dest = 'min_ngal',
                    type = int, default = 10,
                    help = 'Minimum number of members. [Default: 10]')
opts = parser.parse_args()

if not opts.link_r:
    parser.error('Must provide link_r value.')
if not opts.link_z:
    parser.error('Must provide link_z value.')
if not opts.n_pieces:
    parser.error('Must provide n_pieces value.')

# RUN FOF

for i in range(opts.n_pieces):
    code = '/Users/Rowen/Documents/Codes/bin/main'
    fname = opts.dir + 'piece_' + str(i).zfill(2)

    sp.call([code, '-i', fname, '--link_r', str(opts.link_r),
             '--link_z', str(opts.link_z), '--z_min', str(opts.limits[0]),
             '--z_max', str(opts.limits[1]), '--fof_mode', opts.fof_mode,
             '--min_ngal', str(opts.min_ngal)])
