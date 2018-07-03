# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 12:51:14 2018

@author: gregz
"""

import matplotlib
matplotlib.use('agg')
import argparse as ap
import numpy as np
import os.path as op

from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.stats import SigmaClip, biweight_midvariance
from distutils.dir_util import mkpath
from fiber_utils import bspline_x0
from input_utils import setup_logging
from photutils import Background2D, SExtractorBackground, detect_sources
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from utils import biweight_location
from wave_utils import get_new_wave
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.gaussian_process.kernels import ConstantKernel
from sklearn.gaussian_process import GaussianProcessRegressor

import matplotlib.pyplot as plt

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
                        default='/work/03946/hetdex/maverick/red1/reductions')

    parser.add_argument("-op", "--outpath",
                        help='''Outpath for adjusted reductions''',
                        type=str,
                        default='/work/03946/hetdex/maverick/red1/reductions')

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


def set_multi_extensions(outpath, multiname, amps, nfibs, images=[], names=[]):
    for i, amp in enumerate(amps):
        fn = multiname + ('_%s.fits' % amp)
        try:
            F = fits.open(fn)
        except:
            continue
        for name, image in zip(names, images):
            F[name].data = np.array(image[(i*nfibs):((i+1)*nfibs)] * 1.,
                                    dtype='float32')
        F.writeto(op.join(outpath, op.basename(fn)), overwrite=True)


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
    spec, wave, twi, trace, ftf = [np.array(np.vstack(j), dtype='float64')
                                   for j in [spec, wave, twi, tr, ftf]]
    x, y = [np.array(np.hstack(j), dtype='float64') for j in [x, y]]
    return x, y, spec, wave, twi, trace, ftf


def rectify(wave, spec, lims, fac=2.5):
    if wave.ndim == 2:
        N, D = wave.shape
        rect_wave = np.linspace(lims[0], lims[1], int(D*fac))
        rect_spec = np.zeros((N, len(rect_wave)))
        for i in np.arange(N):
            dw = np.diff(wave[i])
            dw = np.hstack([dw[0], dw])
            sel = spec[i] > 1e-3
            if sel.sum() > 10:
                I = interp1d(wave[i][sel], (spec[i] / dw)[sel], kind='quadratic',
                             bounds_error=False, fill_value=-999.)
                rect_spec[i, :] = I(rect_wave)
            else:
                rect_spec[i, :] = 0.0
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


def define_outliers(y, model, neigh=2):
    thresh = 5. * np.sqrt(biweight_midvariance(y - model))
    sel = np.where(np.abs(y - model) > thresh)[0]
    for j in np.arange(1, neigh+1):
        sel = np.union1d(sel, sel + 1)
        sel = np.union1d(sel, sel - 1)
    sel = sel[np.where((sel >= 0) + (sel < len(y)))[0]]
    outlier = np.zeros(y.shape)
    outlier[sel] = True
    return np.array(outlier, dtype=bool)


def fit_bspline(rect_wave, avg, knots=1032):
    B, c = bspline_x0(rect_wave, nknots=knots)
    smooth = np.dot(c, np.linalg.lstsq(c[~avg.mask, :], avg[~avg.mask])[0])
    return np.ma.array(smooth, mask=avg.mask)


def get_avg_spec(wave, spec, twi, lims, mask=None):
    if mask is None:
        mask = np.ones((wave.shape[0],), dtype=bool)
    rect_wave, rect_spec = rectify(wave, spec, lims)
    rect_wave, rect_twi = rectify(wave, twi, lims)
    y = np.ma.array(rect_spec, mask=((rect_spec == 0.) + (rect_spec == -999.)))
    t = np.ma.array(rect_twi, mask=((rect_twi == 0.) + (rect_twi == -999.)))
    fac = np.ma.median(t, axis=1)[:, np.newaxis] / np.ma.median(t)
    norm = y / fac
    avg = biweight_location(norm[mask], axis=(0,))
    smooth = fit_bspline(rect_wave, avg, knots=wave.shape[1])
    return rect_wave, rect_spec, y, norm, avg, smooth, fac


def fit_continuum(wv, sky, sncut=3., skip=1, fil_len=95, func=np.array):
    skym_s = 1. * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    mask = np.zeros(sky.shape, dtype=bool)
    for i in np.arange(5):
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
        I = interp1d(wv[good], sky_sm[good], kind='linear', bounds_error=False,
                     fill_value="extrapolate")
        skym_s[sel] = I(wv[sel])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    mask[sel] = True
    return sky_sm, mask


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


def build_ftf(y, rect_wave, fac, bad, good, nfibs, knots=100):
    mask = np.zeros((y.shape[0],))
    mask[bad] = 1.
    mask = np.array(mask, dtype=bool)
    ftf = []
    for i in np.arange(4):
        xl = i * nfibs
        xh = (i + 1) * nfibs
        X = y[xl:xh] * 1.
        X.mask[mask[xl:xh]] = True
        avg_offset_spec = np.ma.median(X, axis=0)
        smoothed = fit_bspline(rect_wave, avg_offset_spec, knots=knots)
        ftf.append(smoothed + fac[xl:xh])
    return np.vstack(ftf)


def simple_flat_field(X, func=np.abs, fil_len=21, sncut=1.5, skip=2):
    x = np.arange(X.shape[0])
    cont = x * 0.
    outlier = x * 0.
    z = np.ma.median(X, axis=1)
    for i in np.arange(4):
        xl = i * 112
        xh = (i + 1) * 112
        cont[xl:xh], outlier[xl:xh] = fit_continuum(x[xl:xh], z[xl:xh],
                                                    fil_len=fil_len, skip=skip,
                                                    func=func, sncut=sncut)
    return z - cont, cont


def make_frame(xloc, yloc, data, scale=1.,
               seeing_fac=1.5):
    seeing = seeing_fac * scale
    a = len(data)
    x = np.arange(xloc.min()-scale,
                  xloc.max()+1*scale, scale)
    y = np.arange(yloc.min()-scale,
                  yloc.max()+1*scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zimage = xgrid * 0.
    d = np.zeros((a,)+xgrid.shape)
    w = np.zeros((a,)+xgrid.shape)
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            d[:, j, i] = np.sqrt((xloc - xgrid[j, i])**2 +
                                 (yloc - ygrid[j, i])**2)
            w[:, j, i] = np.exp(-1./2.*(d[:, j, i]/seeing)**2)

    sel = np.where((np.abs(data) > 1e-5) * np.isfinite(data))[0]
    ws = w[sel, :, :].sum(axis=0)
    zimage = ((data[sel, np.newaxis, np.newaxis] * w[sel]).sum(axis=0) /
              ws * 1.9)
    return xgrid, ygrid, zimage


def get_sex_background(image, filt_size=21, cols=25):
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(image, (filt_size, cols), filter_size=(1, 1),
                       bkg_estimator=bkg_estimator, sigma_clip=sigma_clip,
                       mask=image.mask, exclude_percentile=100)
    return bkg.background


def setup_GP():
    kernel = (ConstantKernel() + Matern(length_scale=2, nu=3/2) +
              WhiteKernel(noise_level=1.))
    G = GaussianProcessRegressor(alpha=1e-10, copy_X_train=True, kernel=kernel,
                                 n_restarts_optimizer=0, normalize_y=False,
                                 optimizer='fmin_l_bfgs_b', random_state=None)
    return G


def fit_GP(wave, spec, mask):
    G = setup_GP()
    G.fit(wave[mask, np.newaxis], spec[mask])
    return G.predict(wave[:, np.newaxis]), G


def low_cont(x, nfibs):
    z = np.ma.median(X, axis=1)
    z = np.array(z)
    model = z * 0.
    for i in np.arange(4):
        xl = i * nfibs
        xh = (i + 1) * nfibs
        model[xl:xh] = np.percentile(z[xl:xh], 15)
    return z - model, model


def smooth_fiber(X, mask, nfibs):
    z = np.ma.median(X, axis=1)
    z.data[mask] = np.nan
    z = np.array(z)
    x = np.arange(len(z))
    model = z * 0.
    for i in np.arange(4):
        xl = i * nfibs
        xh = (i + 1) * nfibs
        sel = np.isfinite(z[xl:xh])
        G = setup_GP()
        G.fit(x[xl:xh][sel, np.newaxis], z[xl:xh][sel])
        model[xl:xh] = G.predict(x[xl:xh, np.newaxis])
    return model


def make_plot(zimage, xgrid, ygrid, xpos, ypos, good_mask, opath):
    fig = plt.figure(figsize=(6, 6))
    plt.imshow(zimage, origin='lower', interpolation='none', vmin=-15,
               vmax=25, cmap=plt.get_cmap('gray_r'),
               extent=[xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()])
    plt.scatter(xpos[good_mask], ypos[good_mask], marker='x', color='g', s=90)
    plt.scatter(xpos[~good_mask], ypos[~good_mask], marker='x', color='r', s=90)
    fig.savefig(op.join(opath, 'image.png'))


def mask_sources(xgrid, ygrid, xpos, ypos, zimage, sncut=2.0):
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(zimage, (51, 51), filter_size=(1, 1),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       exclude_percentile=100)
    threshold = bkg.background + (sncut * bkg.background_rms)
    kernel = Gaussian2DKernel(2, x_size=5, y_size=5)
    kernel.normalize()
    segm = detect_sources(zimage, threshold, npixels=8, filter_kernel=kernel)
    dist = np.sqrt((xgrid - xpos[:, np.newaxis, np.newaxis])**2 +
                   (ygrid - ypos[:, np.newaxis, np.newaxis])**2)
    fiberloc = np.argmin(dist, axis=0)
    return np.unique(fiberloc[segm.array > 0])

# ['-m', 'multi_317_022_039', '-d', '20180624', '-o', '8', '-e', '1']
args = setup_my_parser(args=None)

if args.instrument == 'virus':
    args.nfibs = 112
if args.instrument == 'lrs2':
    args.nfibs = 140

wave_list = [[3550., 10.], [3735., 10.], [3831., 5.], [3911., 8.],
             [4358, 5.], [4862., 5.], [5085., 5.], [5199., 5.], [5460., 5.]]

for multi in args.multiname:
    args.log.info('Grabbing info for %s' % multi)
    multipath = build_filename(args.rootdir, args.date, args.instrument,
                               args.observation, args.exposure_number, multi)
    xpos, ypos, spec, wave, twi, trace, ftf = get_multi_extensions(multipath,
                                                                   args.amps)
    outpath = build_filename(args.outpath, args.date, args.instrument,
                             args.observation, args.exposure_number, multi)
    outpath = op.dirname(outpath)
    mkpath(outpath)
    args.log.info('Getting average specrum for %s' % multi)
    returned_list = get_avg_spec(wave, spec, twi, args.lims)
    rect_wave, rect_spec, y, norm, avg, smooth, fac = returned_list
    returned_list = get_avg_spec(wave, twi, twi, args.lims)
    rect_wave, rect_twi, y_twi, norm_twi, avg_twi, smooth_twi, fac = returned_list
    if args.recalculate_wavelength:
        newwave = wave * 0.
        for i in np.arange(4):
            xl = i * args.nfibs
            xh = (i+1) * args.nfibs
            args.log.info('Working on the wavelength for fibers: %03d - %03d' %
                          (xl, xh))
            newwave[xl:xh] = get_new_wave(wave[xl:xh], trace[xl:xh],
                                          twi[xl:xh], rect_wave, avg, smooth)
        wave0 = wave * 1.
        wave = newwave * 1.
        returned_list = get_avg_spec(wave, spec, twi, args.lims)
        rect_wave, rect_spec, y, norm, avg, smooth, fac = returned_list
        returned_list = get_avg_spec(wave, twi, twi, args.lims)
        rect_wave, rect_twi, y_twi, norm_twi, avg_twi, smooth_twi, fac = returned_list
    X = rect_spec / rect_twi * smooth_twi / smooth
    flat_field, cont = simple_flat_field(X)
    xgrid, ygrid, zimage = make_frame(xpos, ypos, flat_field)
    mask = mask_sources(xgrid, ygrid, xpos, ypos, zimage)
    good = np.setdiff1d(np.arange(X.shape[0], dtype=int), mask)
    good_mask = np.zeros((X.shape[0],))
    good_mask[good] = 1.
    good_mask = np.array(good_mask, dtype=bool)
    make_plot(zimage * np.ma.median(smooth), xgrid, ygrid, xpos, ypos,
              good_mask, outpath)
    args.log.info('Building fiber to fiber for %s' % multi)
    Y = X * 1.
    faci = []
    SM = X * 1.
    returned_list = get_avg_spec(wave, spec, twi,
                                 args.lims, mask=good_mask)
    rect_wave, rs, y, norm, avg, smooth, fac = returned_list
    Y = (rs - smooth * fac) / smooth
    SM = smooth
    model = smooth_fiber(Y, mask, args.nfibs)[:, np.newaxis]
    Z = (rect_spec - SM * (fac + model)) / SM
    Z.mask[mask] = True
    Z.mask[np.ma.abs(Z) > 0.25] = True
    ftf = rect_spec * 0.
    for i in np.arange(4):
        xl = i * args.nfibs
        xh = (i + 1) * args.nfibs
        back = get_sex_background(Z[xl:xh], 11, 121)
        ftf[xl:xh] = fac[xl:xh] + model[xl:xh] + back
    ZZ = np.ma.median((rect_spec - SM * ftf), axis=1)
    xgrid, ygrid, zimage = make_frame(xpos, ypos, ZZ)
    mask = mask_sources(xgrid, ygrid, xpos, ypos, zimage)
    good = np.setdiff1d(np.arange(X.shape[0], dtype=int), mask)
    good_mask = np.zeros((X.shape[0],))
    good_mask[good] = 1.
    good_mask = np.array(good_mask, dtype=bool)
    make_plot(zimage, xgrid, ygrid, xpos, ypos, good_mask, outpath)
    ftf_new = wave * 0.
    skysub_new = wave * 0.
    sky_new = wave * 0.
    for i in np.arange(wave.shape[0]):
        dw = np.diff(wave[i])
        dw = np.hstack([dw[0], dw])
        I = interp1d(rect_wave, ftf[i], kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        ftf_new[i] = I(wave[i]) * 1.
        J = interp1d(rect_wave, ftf[i]*SM[i], kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        skysub_new[i] = spec[i] - J(wave[i]) * dw
        sky_new[i] = J(wave[i]) * dw

    
    set_multi_extensions(outpath, multipath, args.amps, args.nfibs,
                         images=[ftf_new, sky_new, skysub_new, wave],
                         names=['fiber_to_fiber', 'sky_spectrum',
                                'sky_subtracted', 'wavelength'])


def main():
    pass

if __name__ == main():
    main()
