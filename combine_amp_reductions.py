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


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-f", "--filename",
                    help='''Filename that contains list of files''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''BL=UV, BR=Orange, RL=Red, RR=Farred''',
                    type=str, default=None)

args = parser.parse_args(args=None)

args.log = setup_logging('combine_amp_reductions')
attrs = ['filename', 'side']
for attr in attrs:
    if getattr(args, attr) is None:
        args.log.error('Need a "--%s" argument.' % attr)
        sys.exit(1)


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


def correct_wave(P):
    S = Sky(P.wave, P.spec, P.dar.rect_wave, P.dar.rect_spec,
            P.trace, P.goodfibers)
    newwave = S.wavelength_from_sky()
    return newwave


def rectify(wave, spec, lims, fac=10):
    N, D = wave.shape
    rect_wave = np.linspace(lims[0], lims[1], D*fac)
    rect_spec = np.zeros((N, len(rect_wave)))
    for i in np.arange(N):
        dw = np.diff(wave[i])
        dw = np.hstack([dw[0], dw])
        I = interp1d(wave[i], spec[i] / dw, kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        rect_spec[i, :] = I(rect_wave)
    return rect_wave, rect_spec


def sky_calc(y, goodfibers):
    inds = np.array_split(goodfibers, 14)
    back = y * 0.
    for ind in inds:
        avg = biweight_location(y[ind], axis=(0,))
        xl = ind.min()
        xh = ind.max()+1
        back[xl:xh] = avg
    return back

if args.side == 'RR':
    lims = [8225, 10565]
if args.side == 'RL':
    lims = [6425, 8460]
if args.side == 'BR':
    lims = [4520, 7010]
if args.side == 'BL':
    lims = [3625, 4670]
R = ReduceLRS2(args.filename, args.side)
if args.side[0] == 'R':
    R.dar.spec = 1. * R.oldspec
    R.dar.rectified_dlam = np.abs(np.diff(R.wave_lims)) / (2064.*1.5)
    R.dar.rectify(minwave=R.wave_lims[0], maxwave=R.wave_lims[1])
    R.dar.wave = correct_wave(R)
    R.wave = R.dar.wave * 1.
F = fits.open('ftf_%s.fits' % args.side)
F[0].data = np.array(F[0].data, dtype='float64')
R.ftf = R.wave * 0.
for i in np.arange(R.wave.shape[0]):
    I = interp1d(F[0].data[0], F[0].data[i+1], kind='quadratic',
                 bounds_error=False, fill_value=-999.)
    R.ftf[i] = I(R.wave[i])
rect_wave, rect_spec = rectify(np.array(R.wave, dtype='float64'),
                               np.array(R.oldspec, dtype='float64') / R.ftf,
                               lims)
y = np.ma.array(rect_spec, mask=((rect_spec == 0.) + (rect_spec == -999.)))
back = sky_calc(y, R.goodfibers)
skysub = y * 0.
sky = y * 0.
for i in np.arange(R.wave.shape[0]):
    dw = np.diff(R.wave[i])
    dw = np.hstack([dw[0], dw])
    I = interp1d(rect_wave, back[i], kind='quadratic',
                 bounds_error=False, fill_value=-999.)
    skysub[i] = R.oldspec[i] - I(R.wave[i]) * dw * R.ftf[i]
    sky[i] = I(R.wave[i]) * dw * R.ftf[i]
R.sky = sky * 1.
R.skysub = skysub * 1.
R.ifupos = np.array([R.ifux, R.ifuy])
R.skypos = np.array([R.ra, R.dec])
R.save(image_list=['image', 'error', 'ifupos', 'skypos', 'wave', 'oldspec',
                   'ftf', 'sky', 'skysub'],
       name_list=['image', 'error', 'ifupos', 'skypos', 'wave', 'oldspec',
                  'ftf', 'sky', 'skysub'])
