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


def make_avg_spec(wave, spec, binsize=35, knots=None):
    ''' Make Average spectrum with biweight binning '''
    sel = spec > 0.0
    wave = wave[sel] * 1.
    spec = spec[sel] * 1.
    ind = np.argsort(wave.ravel())
    if wave.ndim == 1:
        N = len(wave)
    else:
        N = wave.shape[0] * wave.shape[1]
    wchunks = np.array_split(wave.ravel()[ind],
                             N / binsize)
    schunks = np.array_split(spec.ravel()[ind],
                             N / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    nspec = np.array([biweight_location(chunk) for chunk in schunks])
    return nwave, nspec


def safe_division(num, denom, eps=1e-8, fillval=0.0):
    good = np.isfinite(denom) * (np.abs(denom) > eps)
    div = num * 0.
    if num.ndim == denom.ndim:
        div[good] = num[good] / denom[good]
        div[~good] = fillval
    else:
        div[:, good] = num[:, good] / denom[good]
        div[:, ~good] = fillval
    return div


def sky_subtract(wave, spec, ftf):
    newspec = safe_division(spec, ftf)
    nwave, nspec = make_avg_spec(np.array(wave, dtype='float64'), 
                                 np.array(newspec, dtype='float64'))
    I = interp1d(nwave, nspec, fill_value='extrapolate', bounds_error=False,
                 kind='quadratic')
    skysub = spec * 0.
    for i in np.arange(wave.shape[0]):
        skysub[i] = spec[i] - I(wave[i]) * ftf[i]
    return skysub


def get_image(fn):
    F = fits.open(fn)
    imagetype = F[0].header['IMAGETYP'].replace(' ', '')
    if imagetype == 'sci':
        S = F['spectrum'].data * 1.
        W = F['wavelength'].data * 1.

#        xarray = np.arange(S.shape[1])
        chunks = np.array_split(S, 20, axis=1)
#        xchunks = np.array([np.mean(x) for x in np.array_split(xarray, 20)])
        avg = np.array([biweight_location(chunk, axis=(1,))
                        for chunk in chunks])
#        I = interp1d(xchunks, avg.swapaxes(0, 1), kind='quadratic',
#                     bounds_error=False, fill_value='extrapolate')
#        norm = I(xarray)
        normavg = biweight_location(avg, axis=(1,))
        divnorm = avg / normavg[:, np.newaxis]
        netnorm = biweight_location(normavg)
        norm = biweight_location(divnorm, axis=(0,)) * netnorm
        return S / norm[:, np.newaxis], W, norm, S
    else:
        return None, None, None, None


def build_residual_frame(dir_list, amp, args, dateb, datee):
    # Create empty lists for the left edge jump, right edge jump, and structure

    sci_list = []
    org_list = []
    norm_list = []
    for directory in dir_list:
        fn = op.join(directory, 'multi_%s_%s.fits' % (args.triplet, amp))
        S, W, N, O = get_image(fn)
        if S is not None:
            sci_list.append(S)
            norm_list.append(N)
            org_list.append(O)

    if not len(sci_list):
        args.log.warning('No reduced frames found for date range given')
        return None
    args.log.info('Number of sci frames from %s-%s for %s: %i' %
                  (dateb, datee, amp, len(sci_list)))
    small_array = np.array(norm_list)
    orig_array = np.array(org_list)
    del org_list
    big_array = np.array(sci_list)
    del sci_list
    func = biweight_location
    mastersci = func(big_array, axis=(0,))

    # Make sky model from average sky
    nwave, nspec = make_avg_spec(np.array(W, dtype='float64'),
                                 np.array(mastersci, dtype='float64'))
    I = interp1d(nwave, nspec, fill_value='extrapolate', bounds_error=False,
                 kind='quadratic')
    ftf = W * 0.
    for fib in np.arange(W.shape[0]):
        ftf[fib] = (mastersci[fib] - I(W[fib])) / I(W[fib])

    # Get average norm
    X = biweight_location(small_array, axis=(1,))[:, np.newaxis]
    norm_of_norms = biweight_location(small_array / X, axis=(0,))
    X = biweight_location(small_array / norm_of_norms[np.newaxis, :],
                          axis=(0,))
    norm_of_norms = biweight_location(small_array / X, axis=(0,))
    master_fiber_to_fiber = ftf + norm_of_norms[:, np.newaxis]
    master_fiber_to_fiber[master_fiber_to_fiber < 0.2] = 0.0

    skysub_list = []
    for i, orig in enumerate(orig_array):
        args.log.info('Making Sky Subtracted frame %i' % (i + 1))
        skysub_list.append(sky_subtract(W, orig, master_fiber_to_fiber))

    sky_array = np.array(skysub_list)
    a, b = master_fiber_to_fiber.shape
    hdu = fits.PrimaryHDU(np.array(master_fiber_to_fiber, dtype='float32'))
    hdu1 = fits.ImageHDU(np.array(W, dtype='float32'))
    hdu2 = fits.ImageHDU(np.array(small_array, dtype='float32'))
    hdu3 = fits.ImageHDU(np.array(sky_array, dtype='float32'))

    mkpath(op.join(args.folder, dateb))
    args.log.info('Writing master_residual_%s_%s.fits' % (args.triplet, amp))
    hdu.header['OBJECT'] = '%s-%s' % (dateb, datee)
    hdu.header['EXTNAME'] = 'fiber_to_fiber'
    hdu1.header['EXTNAME'] = 'wavelength'
    hdu2.header['EXTNAME'] = 'normalization'
    hdu3.header['EXTNAME'] = 'skysub'
    hdulist = fits.HDUList([hdu, hdu1, hdu2, hdu3])
    write_fits(hdulist, op.join(args.folder, dateb,
               'master_residual_%s_%s.fits' % (args.triplet, amp)))

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
