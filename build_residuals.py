# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:32:27 2018

@author: gregz
"""

import glob
import os.path as op
import numpy as np
import sys

from astropy.io import fits
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from utils import biweight_location
from scipy.interpolate import interp1d


def write_fits(hdu, name):
    try:
        hdu.writeto(name, overwrite=True)
    except:
        hdu.writeto(name, clobber=True)


def build_filenames(date, args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    basedir = op.join(args.rootdir, date, args.instrument,
                      args.instrument + '0000*', 'exp*', args.instrument)
    filenames = sorted(glob.glob(op.join(basedir, 'm*_%s_LL.fits' %
                       args.triplet)))
    dirnames = [op.dirname(fn) for fn in filenames]
    unique_dirnames, ind = np.unique(dirnames, return_index=True)
    return list(unique_dirnames)


def get_image(fn):
    F = fits.open(fn)
    imagetype = F[0].header['IMAGETYP'].replace(' ', '')
    if imagetype == 'sci':
        S = F['spectrum'].data * 1.
        W = F['spectrum'].data * 1.

        xarray = np.arange(S.shape[1])
        chunks = np.array_split(S, 20, axis=1)
        xchunks = np.array([np.mean(x) for x in np.array_split(xarray, 20)])
        avg = np.array([biweight_location(chunk, axis=(1,))
                        for chunk in chunks])
        I = interp1d(xchunks, avg.swapaxes(0, 1), kind='quadratic',
                     bounds_error=False, fill_value='extrapolate')
        norm = I(xarray)
        return S / norm, W
    else:
        return None, None


def build_residual_frame(dir_list, amp, args, dateb, datee):
    # Create empty lists for the left edge jump, right edge jump, and structure

    sci_list = []
    for directory in dir_list:
        fn = op.join(directory, 'multi_%s_%s.fits' % (args.triplet, amp))
        S, W = get_image(fn)
        if S is not None:
            sci_list.append(S)

    if not len(sci_list):
        args.log.warning('No reduced frames found for date range given')
        return None
    args.log.info('Number of sci frames from %s-%s for %s: %i' % (date_begin,
                                                                  date_end,
                                                                  amp,
                                                                  len(sci_list)))
    big_array = np.array(sci_list)
    func = biweight_location
    mastersci = func(big_array, axis=(0,))

    a, b = mastersci.shape
    hdu = fits.PrimaryHDU(np.array(mastersci, dtype='float32'))
    hdu1 = fits.ImageHDU(np.array(W, dtype='float32'))
    mkpath(op.join(args.folder, dateb))
    args.log.info('Writing master_residual_%s_%s.fits' % (args.triplet, amp))
    hdu.header['OBJECT'] = '%s-%s' % (dateb, datee)
    hdu.header['EXTNAME'] = 'spectrum'
    hdu1.header['EXTNAME'] = 'wavelength'
    hdulist = fits.HDUList([hdu, hdu1])
    write_fits(hdulist, op.join(args.folder, dateb, 'master_residual_%s_%s.fits' %
               (args.triplet, amp)))

parser = setup_parser()
parser.add_argument("-f", "--folder",
                    help='''Output folder''',
                    type=str, default='residuals')

parser.add_argument("-m", "--maxnum",
                    help='''Maximum number of bias frames in masterbias''',
                    type=int, default=100)

parser.add_argument("-tr", "--triplet",
                    help='''Triplet of the specid, ifuslot, ifuid''',
                    type=str, default=None)

amps = ['LL', 'LU', 'RL', 'RU']

args = parser.parse_args(args=None)
args.log = setup_logging(logname='build_master_bias')
args = set_daterange(args)
if args.triplet is None:
    args.log.error('Please set the "--triplet" argument.')
    sys.exit(1)

filenames = []
for date in args.daterange:
    date = '%04d%02d%02d' % (date.year, date.month, date.day)
    filenames = filenames + build_filenames(date, args)

date_begin = args.daterange[0]
date_end = args.daterange[-1]
date_begin = '%04d%02d%02d' % (date_begin.year, date_begin.month,
                               date_begin.day)
date_end = '%04d%02d%02d' % (date_end.year, date_end.month, date_end.day)

args.log.info('Length of filenames found for %s-%s: %i' % (date_begin, 
                                                           date_end,
                                                           len(filenames)))
if (len(filenames) % args.maxnum) == 0:
    nbins = len(filenames) / args.maxnum
else:
    nbins = len(filenames) / args.maxnum + 1
if nbins == 0:
    args.log.warning('No files found for %s on %s-%s' % (args.triplet,
                                                         date_begin,
                                                         date_end))
    sys.exit(1)

chunks = np.array_split(filenames, nbins)
for chunk in chunks:
    bases = [chunk[0], chunk[-1]]
    dates = []
    for base in bases:
        base0 = str(base)
        for i in np.arange(4):
            base0 = op.dirname(base0)
        dates.append(op.basename(base0))
    for amp in amps:
        build_residual_frame(chunk, amp, args, dates[0], dates[1])
