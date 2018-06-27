# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 12:51:14 2018

@author: gregz
"""

import argparse as ap
import numpy as np
import os.path as op

from astropy.io import fits
from astropy.stats import SigmaClip, biweight_midvariance
from fiber_utils import bspline_x0
from input_utils import setup_logging
from photutils import Background2D, BiweightLocationBackground
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from utils import biweight_location
from wave_utils import get_new_wave


def setup_my_parser(args=None):
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-e", "--exposure_number",
                        help='''Exposure number, 10''',
                        type=int, default=None)

    parser.add_argument("-m", "--multiname",
                        help='''multi* base name leaving off "_LL.fits"''',
                        type=str, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str,
                        default='work/03946/hetdex/maverick/red1/reductions')

    parser.add_argument("-in", "--instrument",
                        help='''Instrument, e.g., virus''',
                        type=str, default='virus')

    parser.add_argument("-a", "--amps",
                        help='''amplifiers to use''',
                        type=str, default='LL, LU, RU, RL')

    parser.add_argument("-wl", "--wave_lims",
                        help='''wavelength limits of the instrume''',
                        type=str, default='3470, 5530')

    parser.add_argument("-rc", "--recalculate_wavelength",
                        help='''recalculate_wavelength''',
                        action="count", default=0)
    args = parser.parse_args(args=args)
    args.log = setup_logging('great_code_ever')

    attr = ['date', 'observation', 'exposure_number', 'multiname']
    for att in attr:
        if getattr(args, att) is None:
            args.log.error('Please set "--%s" argument.' % att)
            return None

    args.multiname = [j.replace(' ', '') for j in args.multiname.split(',')]
    args.amps = [j.replace(' ', '') for j in args.amps.split(',')]
    args.lims = [float(j.replace(' ', '')) for j in args.wave_lims.split(',')]
    return args


def build_filename(rootdir, date, instrument, obs, expn, multiname):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    filename = op.join(rootdir, date, instrument,
                       '%s%07d' % (instrument, int(obs)),
                       'exp%02d' % expn, instrument,
                       multiname)
    return filename


def get_multi_extensions(multiname, amps):
    x, y, spec, wave, twi, tr, ftf = ([], [], [], [], [], [], [])
    for amp in amps:
        fn = multiname + ('_%s.fits' % amp)
        try:
            F = fits.open(fn)
        except:
            continue
        try:
            x.append(F['ifupos'].data[:, 0])
            y.append(F['ifupos'].data[:, 1])
        except:
            x.append(np.ones((112,)))
            y.append(np.ones((112,)))
        spec.append(F['spectrum'].data)
        try:
            twi.append(F['twi_spectrum'].data)
        except:
            twi.append(F['spectrum'].data)
        wave.append(F['wavelength'].data)
        ftf.append(F['fiber_to_fiber'].data)
        if amp in ['LU', 'RL']:
            addtr = F[0].data.shape[0]
        else:
            addtr = 0.
        tr.append(F['trace'].data + addtr)
    spec, wave, twi, trace, ftf = [np.vstack(j)
                                   for j in [spec, wave, twi, tr, ftf]]
    x, y = [np.hstack(j) for j in [x, y]]
    return x, y, spec, wave, twi, trace, ftf


def rectify(wave, spec, lims, fac=2.5):
    if wave.ndim == 2:
        N, D = wave.shape
        rect_wave = np.linspace(lims[0], lims[1], int(D*fac))
        rect_spec = np.zeros((N, len(rect_wave)))
        for i in np.arange(N):
            dw = np.diff(wave[i])
            dw = np.hstack([dw[0], dw])
            I = interp1d(wave[i], spec[i] / dw, kind='quadratic',
                         bounds_error=False, fill_value=-999.)
            rect_spec[i, :] = I(rect_wave)
    else:
        D = len(wave)
        rect_wave = np.linspace(lims[0], lims[1], int(D*fac))
        rect_spec = np.zeros((len(rect_wave, )))
        dw = np.diff(wave)
        dw = np.hstack([dw[0], dw])
        I = interp1d(wave, spec / dw, kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        rect_spec = I(rect_wave)
    return rect_wave, rect_spec


def fit_bspline(rect_wave, avg):
    B, c = bspline_x0(rect_wave, nknots=wave.shape[1])
    smooth = np.dot(c, np.linalg.lstsq(c[~avg.mask, :], avg[~avg.mask])[0])
    return np.ma.array(smooth, mask=avg.mask)


def get_avg_spec(wave, spec, twi, lims):
    rect_wave, rect_spec = rectify(wave, spec, lims)
    rect_wave, rect_twi = rectify(wave, twi, lims)
    y = np.ma.array(rect_spec, mask=((rect_spec == 0.) + (rect_spec == -999.)))
    t = np.ma.array(rect_twi, mask=((rect_twi == 0.) + (rect_twi == -999.)))
    fac = np.ma.median(t, axis=1)[:, np.newaxis] / np.ma.median(t)
    norm = y / fac
    avg = biweight_location(norm, axis=(0,))
    smooth = fit_bspline(rect_wave, avg)
    return rect_wave, rect_spec, y, norm, avg, smooth, fac


def fit_continuum(wv, sky, sncut=3., skip=1, fil_len=95, func=np.array):
    skym_s = 1. * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    for i in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        mad = np.sqrt(biweight_midvariance(sky-sky_sm))
        outlier = func(sky - sky_sm) > sncut * mad
        sel = np.where(outlier)[0]
        for j in np.arange(1, skip+1):
            sel = np.union1d(sel, sel + 1)
            sel = np.union1d(sel, sel - 1)
        sel = np.sort(np.unique(sel))
        sel = sel[skip:-skip]
        good = np.setdiff1d(allind, sel)
        skym_s = 1.*sky
        skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    return sky_sm, sel


def longest_consecutive_mask(outliers):
    tot = 0
    cnt = 0
    for i in np.arange(1, len(outliers)):
        if outliers[i] == (outliers[i-1] + 1):
            cnt += 1
        else:
            tot = np.max([tot, cnt])
            cnt = 0
    return tot


def build_ftf(norm, smooth, fac):
    y = (norm - smooth) / smooth
    y[y.mask] = 0.0
    y.mask[y.data == 0.0] = True
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = BiweightLocationBackground()
    bkg1 = Background2D(y.data[:224], (15, 251), filter_size=(1, 1),
                        sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    bkg2 = Background2D(y.data[224:], (15, 251), filter_size=(1, 1),
                        sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    background = np.vstack([bkg1.background, bkg2.background])
    z = np.ma.median(y - background, axis=1)
    x = np.arange(len(z))
    cont, outliers = fit_continuum(x, z)
    L = 2 * longest_consecutive_mask(outliers)
    if not (L & 1):
        L += 1
    L = np.max([L, 15])
    mask = np.zeros(y.data.shape, dtype=bool)
    mask[outliers] = True
    mask[y.data == 0.0] = True
    bkg1 = Background2D(y.data[:224], (L, 51), filter_size=(1, 1),
                        mask=mask[:224], sigma_clip=sigma_clip,
                        bkg_estimator=bkg_estimator)
    bkg2 = Background2D(y.data[224:], (L, 51), filter_size=(1, 1),
                        mask=mask[224:], sigma_clip=sigma_clip,
                        bkg_estimator=bkg_estimator)
    background = np.vstack([bkg1.background, bkg2.background])
    return fac + background, x, z, cont, outliers, bkg1


args = setup_my_parser(args=['-d', '20180624', '-o', '8', '-e', '1', '-m',
                             'multi_321_026_035'])

for multi in args.multiname:
    multipath = build_filename(args.rootdir, args.date, args.instrument,
                               args.observation, args.exposure_number, multi)
    xpos, ypos, spec, wave, twi, trace, ftf = get_multi_extensions(multipath,
                                                                   args.amps)
    returned_list = get_avg_spec(wave, spec, twi, args.lims)
    rect_wave, rect_spec, y, norm, avg, smooth, fac = returned_list
    if args.recalculate_wavelength:
        newwave = wave * 0.
        args.log.info('Working on the wavelength for side L')
        newwave[:224] = get_new_wave(wave[:224], trace[:224], twi[:224],
                                     rect_wave, avg, smooth)
        args.log.info('Working on the wavelength for side R')
        newwave[224:] = get_new_wave(wave[224:], trace[224:], twi[224:],
                                     rect_wave, avg, smooth)
        returned_list = get_avg_spec(newwave, spec, twi, args.lims)
        rect_wave, rect_spec, y, norm, avg, smooth, fac = returned_list
    ftf, x, z, cont, outliers, bkg = build_ftf(norm, smooth, fac)
    new_norm = y / ftf
    new_avg = biweight_location(new_norm, axis=(0,))
    new_smooth = fit_bspline(rect_wave, new_avg)
    new_ftf, x, z, cont, outliers, bkg = build_ftf(new_norm, new_smooth, ftf)
    fits.PrimaryHDU(new_ftf.data).writeto('test.fits', overwrite=True)
    fits.PrimaryHDU((y - new_ftf*new_smooth).data).writeto('test_res.fits',
                                                           overwrite=True)


def main():
    pass

if __name__ == main():
    main()
