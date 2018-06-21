# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:02:02 2018

@author: gregz
"""

import numpy as np
import sys
from reducelrs2 import ReduceLRS2
from astropy.io import fits
from scipy.interpolate import interp1d
from input_utils import setup_logging
import argparse as ap
from skysubtraction import Sky
from utils import biweight_location
import glob
import os.path as op


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-s", "--side",
                    help='''UV=BL, Orange=BR, Red=RL, Farred=RR''',
                    type=str, default=None)
parser.add_argument("-o", "--outname",
                    help='''Name of the fiber to fiber frame output''',
                    type=str, default='ftf.fits')
parser.add_argument("-d", "--daterange",
                    help='''Inclusive date range for building fib to fib''',
                    type=str, default='20180101, 20180501')
parser.add_argument("-e", "--exposuretimerange",
                    help='''Exposure range for building fib to fib''',
                    type=str, default='120, 7200')
parser.add_argument("-r", "--reductiondir",
                    help='''Name of the reduction directory''',
                    type=str, default='reductions')
parser.add_argument("-n", "--nbins",
                    help='''Number of bins for smoothing''',
                    type=int, default=15)

args = parser.parse_args(args=None)

args.log = setup_logging('new_fiber_to_fiber')
attrs = ['side']
for attr in attrs:
    if getattr(args, attr) is None:
        args.log.error('Need a "--%s" argument.' % attr)
        sys.exit(1)

args.dates = [x.replace(' ', '') for x in args.daterange.split(',')]
args.exptimes = [x.replace(' ', '') for x in args.exposuretimerange.split(',')]
args.exptimes = [float(x) for x in args.exptimes]
side_dict = {'BL': ['056', ['LL', 'LU'], [3625., 4670.], '066'],
             'BR': ['056', ['RL', 'RU'], [4520., 7010.], '066'],
             'RL': ['066', ['LL', 'LU'], [6425., 8460.], '056'],
             'RR': ['066', ['RL', 'RU'], [8225., 10565.], '056']}
searchname = op.join(args.reductiondir, '*/lrs2/*/*/lrs2/m*%s*.fits' %
                                        (side_dict[args.side][0]))
filenames = sorted(glob.glob(searchname))
filelist = []
for fn in filenames:
    date = fn.split('/')[-6]
    if (date >= args.dates[0]) and (date <= args.dates[1]):
        cond1 = side_dict[args.side][1][0] in fn
        cond2 = side_dict[args.side][1][1] in fn
        if cond1 or cond2:
            F = fits.open(fn)
            cond1 = F[0].header['EXPTIME'] >= args.exptimes[0]
            cond2 = F[0].header['EXPTIME'] <= args.exptimes[1]
            cond3 = side_dict[args.side][3] in F[0].header['OBJECT']
            if cond1 and cond2 and cond3:
                filelist.append(fn[:-8])

filelist = np.unique(filelist)
print(len(filelist))


def smooth(rect_wave, spec_array, n=15):
    chunks = np.array_split(spec_array, n, axis=1)
    mvals = [np.ma.median(chunk, axis=1) for chunk in chunks]
    waves = np.array_split(rect_wave, n)
    mwave = [np.median(wave) for wave in waves]
    ftf = spec_array * 0.
    for i in np.arange(spec_array.shape[0]):
        I = interp1d(mwave, [mval[i] for mval in mvals], kind='quadratic',
                     bounds_error=False, fill_value="extrapolate")
        ftf[i, :] = I(rect_wave)
    return ftf


def build_ftf(rect_wave, spec_array, n=15):
    chunks = np.array_split(spec_array, n, axis=2)
    mvals = [np.ma.median(chunk, axis=(0, 2)) for chunk in chunks]
    waves = np.array_split(rect_wave, n)
    mwave = [np.median(wave) for wave in waves]
    ftf = spec_array[0] * 0.
    for i in np.arange(spec_array.shape[1]):
        I = interp1d(mwave, [mval[i] for mval in mvals], kind='quadratic',
                     bounds_error=False, fill_value="extrapolate")
        ftf[i, :] = I(rect_wave)
    return ftf


def correct_wave(P):
    S = Sky(P.wave, P.spec, P.dar.rect_wave, P.dar.rect_spec,
            P.trace, P.goodfibers)
    newwave = S.wavelength_from_sky()
    return newwave


def rectify(wave, spec, lims, fac=2.5):
    N, D = wave.shape
    rect_wave = np.linspace(lims[0], lims[1], int(D*fac))
    rect_spec = np.zeros((N, len(rect_wave)))
    for i in np.arange(N):
        dw = np.diff(wave[i])
        dw = np.hstack([dw[0], dw])
        I = interp1d(wave[i], spec[i] / dw, kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        rect_spec[i, :] = I(rect_wave)
    return rect_wave, rect_spec


spec_list = []
base_list = []
wave_list = []

for filebase in filelist:
    R = ReduceLRS2(filebase, args.side)
#    if args.side[0] == 'R':
#        R.dar.spec = 1. * R.oldspec
#        R.dar.rectified_dlam = np.abs(np.diff(R.wave_lims)) / (2064.*1.5)
#        R.dar.rectify(minwave=R.wave_lims[0], maxwave=R.wave_lims[1])
#        R.dar.wave = correct_wave(R)
    rect_wave, rect_spec = rectify(np.array(R.wave, dtype='float64'),
                                   np.array(R.oldspec, dtype='float64'),
                                   side_dict[args.side][2])
    y = np.ma.array(rect_spec, mask=((rect_spec == 0.) + (rect_spec == -999.)))
    norm = (y / biweight_location(y, axis=(1,))[:, np.newaxis] *
            biweight_location(y))
    avg = biweight_location(norm, axis=(0,))
    spec_list.append(y / avg)
    wave_list.append(R.wave * 1.)

y = np.array(spec_list)
y = np.ma.array(y, mask=(np.isnan(y) + (y <= 0)))
ftf = build_ftf(rect_wave, y, n=args.nbins)
FtF = np.vstack([rect_wave, ftf])
F = fits.PrimaryHDU(np.array(FtF, dtype='float32')).writeto(args.outname,
                                                            overwrite=True)
