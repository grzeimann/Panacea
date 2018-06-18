# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:02:02 2018

@author: gregz
"""

import numpy as np
import sys
import os.path as op
from reducelrs2 import ReduceLRS2
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d
from input_utils import setup_logging
import argparse as ap
from skysubtraction import Sky
from utils import biweight_location, biweight_midvariance
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.modeling.models import Moffat2D, Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-f", "--filename",
                    help='''Filename that contains list of files''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''BL=UV, BR=Orange, RL=Red, RR=Farred''',
                    type=str, default=None)
parser.add_argument("-g", "--fibergroup",
                    help=''' Size of the fiber group for sky subtraction''',
                    type=int, default=20)

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


def safe_division(num, denom, eps=1e-8, fillval=0.0):
    good = np.isfinite(denom) * (np.abs(denom) > eps)
    div = num * 0.
    div[good] = num[good] / denom[good]
    div[~good] = fillval
    return div


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


def sky_calc(y, goodfibers, nbins=14):
    inds = np.array_split(goodfibers, nbins)
    back = y * 0.
    for ind in inds:
        avg = biweight_location(y[ind], axis=(0,))
        xl = ind.min()
        xh = ind.max()+1
        back[xl:xh] = avg
    return back


def gather_sn_fibers(fibconv, noise, cols):
    hightolow = np.argsort(np.median(fibconv[:, cols], axis=1))[::-1]
    s = 0.
    ss = np.zeros((len(cols),))
    nn = noise[cols]
    inds = []
    for ind in hightolow:
        news = fibconv[ind, cols] + ss
        newn = np.sqrt(nn**2 + noise[cols]**2)
        rat = np.median(news / newn)
        if rat > s:
            nn = newn
            ss = news
            s = rat
            inds.append(ind)
        else:
            continue
    return inds, s


def find_centroid(image, x, y):
    G = Gaussian2D()
    fit = LevMarLSQFitter()(G, x, y, image)
    return fit.x_mean.value, fit.y_mean.value


def build_big_fiber_array(P):
    u = np.unique(P.ifuy)
    NX = []
    NY = []
    for ui in u:
        X = np.sort(P.ifux[P.ifuy == ui])
        lx = X[-1] - X[0]
        dx = X[1] - X[0]
        NX.append(np.hstack([X - lx - dx, X, X + lx + dx]))
        NY.append(ui * np.ones(NX[-1].shape))
    NX.append(NX[-2])
    NY.append(NY[-2]+0.59*2.)
    NX = np.hstack(NX)
    NY = np.hstack(NY)
    ly = NY[-1] - NY[0]
    dy = u[1] - u[0]
    NX = np.hstack([NX, NX, NX])
    NY = np.hstack([NY - ly - dy, NY, NY + ly + dy])
    return NX, NY


def flux_correction(wave, loc, P, inds, dar_table, rect_spec,
                    rect_sky, alpha=3.5, gamma=1.5):
    X = interp1d(dar_table['wave'], dar_table['x_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    Y = interp1d(dar_table['wave'], dar_table['y_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    frac = wave * 0.
    SF = wave * 0.
    SS = wave * 0.
    NX, NY = build_big_fiber_array(P)
    PSF = Moffat2D(amplitude=1., x_0=0., y_0=0., alpha=alpha, gamma=gamma)
    for i in np.arange(len(wave)):
        x = loc[1] + X(wave[i]) - X(loc[0])
        y = loc[2] + Y(wave[i]) - Y(loc[0])
        PSF.x_0.value = x
        PSF.y_0.value = y
        total = PSF(NX, NY).sum()
        wei = PSF(P.ifux[inds], P.ifuy[inds])
        SF[i] = np.sum(rect_spec[inds, i] * wei) * total / wei.sum()
        SS[i] = np.sum(rect_sky[inds, i] * wei) * total / wei.sum()
        frac[i] = wei.sum() / total
    return frac, SF, SS


if args.side == 'RR':
    lims = [8225, 10565]
if args.side == 'RL':
    lims = [6425, 8460]
if args.side == 'BR':
    lims = [4520, 7010]
if args.side == 'BL':
    lims = [3625, 4670]


def main():
    R = ReduceLRS2(args.filename, args.side)
    R.get_mirror_illumination()
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
                                   np.array(R.oldspec, dtype='float64') /
                                   R.ftf, lims, fac=2.5)

    y = np.ma.array(rect_spec, mask=((rect_spec == 0.) + (rect_spec == -999.)))

    for i in np.arange(2):
        back = sky_calc(y, R.goodfibers, nbins=(R.wave.shape[0] /
                                                args.fibergroup))
        G = Gaussian1DKernel(1.5)
        fibconv = rect_spec * 0.
        for i in np.arange(R.wave.shape[0]):
            fibconv[i] = convolve(rect_spec[i] - back[i], G)
        noise = biweight_midvariance(fibconv, axis=(0,))
        R.signoise = fibconv / noise
        S = np.nanmedian(R.signoise, axis=1)
        N = biweight_midvariance(S)
        R.goodfibers = np.where((S/N) < 3.)[0]

    skysub = R.wave * 0.
    sky = R.wave * 0.
    for i in np.arange(R.wave.shape[0]):
        dw = np.diff(R.wave[i])
        dw = np.hstack([dw[0], dw])
        I = interp1d(rect_wave, back[i], kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        skysub[i] = R.oldspec[i] - I(R.wave[i]) * dw * R.ftf[i]
        sky[i] = I(R.wave[i]) * dw * R.ftf[i]

    R.sky = sky * 1.
    R.skysub = safe_division(skysub, R.ftf)
    R.ifupos = np.array([R.ifux, R.ifuy]).swapaxes(0, 1)
    R.skypos = np.array([R.ra, R.dec]).swapaxes(0, 1)
    R.skynorm = safe_division(R.sky, R.ftf)

    T = Table.read('response_%s.dat' % args.side,
                   format='ascii.fixed_width_two_line')
    I = interp1d(T['wave'], T['response'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    R.flam = R.wave * 0.
    R.slam = R.wave * 0.
    for i in np.arange(R.wave.shape[0]):
        response = I(R.wave[i])
        R.flam[i] = R.skysub[i] / R.exptime / R.area * response
        R.slam[i] = R.skynorm[i] / R.exptime / R.area * response

    rect_wave, rect_spec = rectify(np.array(R.wave, dtype='float64'),
                                   np.array(R.flam, dtype='float64'),
                                   lims, fac=1.0)
    rect_wave, rect_sky = rectify(np.array(R.wave, dtype='float64'),
                                  np.array(R.slam, dtype='float64'),
                                  lims, fac=1.0)

    R.define_good_fibers()
    G = Gaussian1DKernel(1.5)
    fibconv = rect_spec * 0.
    for i in np.arange(R.wave.shape[0]):
        fibconv[i] = convolve(rect_spec[i], G)
    noise = biweight_midvariance(fibconv, axis=(0,))
    inds = np.array_split(np.arange(len(rect_wave)), 20)
    v = np.argmax([np.sum(np.nanmedian(chunk, axis=1))
                   for chunk in np.array_split(fibconv / noise, 20, axis=1)])
    sn_image = np.nanmedian(np.array_split(fibconv / noise, 20, axis=1)[v],
                            axis=1)
    wv = np.median(np.array_split(rect_wave, 20)[v])
    xc, yc = find_centroid(sn_image, R.ifux, R.ifuy)
    print(v, xc, yc)

    fibinds, s = gather_sn_fibers(fibconv, noise, inds[v])
    dar_table = Table.read('dar_%s.dat' % args.side,
                           format='ascii.fixed_width_two_line')

    frac, R.flux, R.skyflux = flux_correction(rect_wave, [wv, xc, yc], R,
                                              fibinds, dar_table,
                                              rect_spec, rect_sky)
    print(len(fibinds), s, np.median(frac))

    # R.flux = rect_spec[np.array(fibinds, dtype=int), :].sum(axis=0) / frac
    # R.skyflux = rect_sky[np.array(fibinds, dtype=int), :].sum(axis=0) / frac
    R.fluxerror = noise * np.sqrt(len(fibinds)) / frac
    R.save(image_list=['image_name', 'error', 'ifupos', 'skypos', 'wave',
                       'oldspec', 'ftf', 'sky', 'skysub'],
           name_list=['image', 'error', 'ifupos', 'skypos', 'wave', 'oldspec',
                      'ftf', 'sky', 'skysub'])
    names = ['wavelength', 'F_lambda', 'e_F_lambda', 'Sky_lambda',
             'e_Sky_lambda']
    hdu = fits.PrimaryHDU(np.array([rect_wave, R.flux, R.fluxerror, R.skyflux],
                                   dtype='float32'))
    hdu.header['waveunit'] = 'A'
    hdu.header['fluxunit'] = 'ergs/s/cm2/A'
    for i, name in enumerate(names):
        hdu.header['ROW%i' % (i+1)] = name
    for key in R.header.keys():
        if key in hdu.header:
            continue
        hdu.header[key] = R.header[key]
    hdu.writeto(op.join(R.path, 'spectrum_%s.fits' % args.side),
                overwrite=True)

if __name__ == '__main__':
    main()
