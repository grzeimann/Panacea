# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:02:02 2018

@author: gregz
"""

import matplotlib
matplotlib.use('agg')

import argparse as ap
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import sys
import cosmics

from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Moffat2D
from astropy.table import Table
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from copy import copy
from fiber_utils import bspline_x0
from input_utils import setup_logging
from photutils import detect_sources
from reducelrs2 import ReduceLRS2
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.gaussian_process.kernels import ConstantKernel
from sklearn.gaussian_process import GaussianProcessRegressor
from utils import biweight_location, biweight_midvariance
from wave_utils import get_new_wave, get_red_wave, get_single_shift


def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()
parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-f", "--filename",
                    help='''Filename that contains list of files''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''blue for LRS2-B and red for LRS2-R''',
                    type=str, default='blue')
parser.add_argument("-rc", "--recalculate_wavelength",
                    help='''recalculate_wavelength''',
                    action="count", default=0)
parser.add_argument("-em", "--emission",
                    help='''Find emission line object?''',
                    action="count", default=0)
parser.add_argument("-es", "--extract_side",
                    help='''blue for LRS2-B and red for LRS2-R''',
                    type=str, default='orange')
parser.add_argument("-we", "--wave_extract",
                    help='''blue for LRS2-B and red for LRS2-R''',
                    type=float, default=None)
args = parser.parse_args(args=None)

args.log = setup_logging('combine_amp_reductions')
attrs = ['filename', 'side']
for attr in attrs:
    if getattr(args, attr) is None:
        args.log.error('Need a "--%s" argument.' % attr)
        sys.exit(1)
args.side = args.side.lower()


def make_avg_spec(wave, spec, binsize=35, knots=None):
    if knots is None:
        knots = wave.shape[1]
    ind = np.argsort(wave.ravel())
    N, D = wave.shape
    wchunks = np.array_split(wave.ravel()[ind],
                             N * D / binsize)
    schunks = np.array_split(spec.ravel()[ind],
                             N * D / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    B, c = bspline_x0(nwave, nknots=knots)
    nspec = np.array([biweight_location(chunk) for chunk in schunks])
    sol = np.linalg.lstsq(c, nspec)[0]
    smooth = np.dot(c, sol)
    nwave, nind = np.unique(nwave, return_index=True)
    return nwave, smooth[nind]


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
        if rat > (s+0.5):
            nn = newn
            ss = news
            s = rat
            inds.append(ind)
        else:
            continue
    return inds, s


def find_centroid(image, x, y, B):
    G = Moffat2D()
    G.alpha.value = 3.5
    G.alpha.fixed = True
    fit = LevMarLSQFitter()(G, x, y, image)
    signal_to_noise = fit.amplitude.value / biweight_midvariance(image)
    d = np.sqrt((x - fit.x_0.value)**2 + (y - fit.y_0.value)**2)
    ratio = fit(x, y) / B
    ind = np.argsort(ratio)
    dthresh = np.interp(.01, ratio[ind], d[ind])
    return (fit.x_0.value, fit.y_0.value, fit.alpha.value, fit.gamma.value,
            fit.fwhm, signal_to_noise, dthresh)


def build_weight_matrix(x, y, sig=1.5):
    d = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    G = np.exp(-0.5 * (d / sig)**2)
    G = G / G.sum(axis=0)[:, np.newaxis]
    return G.swapaxes(0,1)


def clean_cosmics(rect_spec):
    G = np.array([-.25, -.25, 1., -.25, -.25]).reshape(5,1)
    S = convolve(rect_spec, G, normalize_kernel=False)
    N = biweight_midvariance(S, axis=(0,))
    mask = 0. * rect_spec
    mask[(S / N) > 5.] = -1.
    print('[COSMICS] Number of cosmics found: %i' % int(-1.*mask.sum()))
    return mask


def mask_skylines_cosmics(wave, rect_spec, name):
    mask1 = rect_spec * 0.
    if op.exists(op.join(DIRNAME, 'lrs2_config', '%s_skylines.dat' % name)):
        T = Table.read(op.join(DIRNAME, 'lrs2_config', '%s_skylines.dat' % name),
                       format='ascii.fixed_width_two_line')
        for w in T['wavelength']:
            mask1[:, np.abs(wave - w) < 6.] = -1.
    mask2 = clean_cosmics(rect_spec)
    mask = (mask1 + mask2) < 0
    return mask


def convolve_spatially(x, y, spec, wave, name, sig_spatial=0.7, sig_wave=1.5):
    W = build_weight_matrix(x, y, sig=sig_spatial)
    mask = mask_skylines_cosmics(wave, spec, name)
    Z = spec * 1.
    Z[mask] = np.nan
    G = Gaussian1DKernel(sig_wave)
    for i in np.arange(spec.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment='fill', fill_value=0.0)
    for i in np.arange(spec.shape[1]):
        Z[:, i] = np.dot(Z[:, i], W)
    return Z, mask


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


def get_x_y_lambda(det_ind, other_ind, detwv, otherwv,
                   det_xc, det_yc, sides):
    side = sides[det_ind]
    dar_table = Table.read(op.join(DIRNAME, 'lrs2_config', 'dar_%s.dat' % side),
                           format='ascii.fixed_width_two_line')
    X = interp1d(dar_table['wave'], dar_table['x_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    Y = interp1d(dar_table['wave'], dar_table['y_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    xoff, yoff = (det_xc - X(detwv), det_yc - Y(detwv))
    side = sides[other_ind]
    dar_table = Table.read(op.join(DIRNAME, 'lrs2_config', 'dar_%s.dat' % side),
                           format='ascii.fixed_width_two_line')
    X = interp1d(dar_table['wave'], dar_table['x_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    Y = interp1d(dar_table['wave'], dar_table['y_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    return xoff + X(otherwv), yoff + Y(otherwv)


def make_plot(zimage, xgrid, ygrid, xpos, ypos, good_mask, opath, side):
    fig = plt.figure(figsize=(6, 6))
    plt.imshow(zimage, origin='lower', interpolation='none',
               norm=ImageNormalize(stretch=AsinhStretch()),
               cmap=plt.get_cmap('gray_r'),
               extent=[xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()])
    plt.scatter(xpos[good_mask], ypos[good_mask], marker='x', color='g', s=90)
    plt.scatter(xpos[~good_mask], ypos[~good_mask], marker='x', color='r',
                s=90)
    plt.axis([xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()])
    fig.savefig(op.join(opath, 'image_%s.png' % side))


def mask_sources(xgrid, ygrid, xpos, ypos, zimage, sncut=2.0):
    threshold = (biweight_location(zimage) +
                 sncut * biweight_midvariance(zimage))
    kernel = Gaussian2DKernel(2, x_size=5, y_size=5)
    kernel.normalize()
    segm = detect_sources(zimage, threshold, npixels=8, filter_kernel=kernel)
    dist = np.sqrt((xgrid - xpos[:, np.newaxis, np.newaxis])**2 +
                   (ygrid - ypos[:, np.newaxis, np.newaxis])**2)
    fiberloc = np.argmin(dist, axis=0)
    return np.unique(fiberloc[segm.array > 0])


def make_frame(xloc, yloc, data, scale=0.25,
               seeing_fac=2.5):
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


def smooth_fiber(X, mask, nfibs, wave_sel=None):
    z = biweight_location(X, axis=(1,))
    x = np.arange(len(z))
    z[mask] = np.nan
    model = z * 0.
    for i in np.arange(2):
        xl = i * nfibs
        xh = (i + 1) * nfibs
        sel = np.isfinite(z[xl:xh])
        G = setup_GP()
        G.fit(x[xl:xh][sel, np.newaxis], z[xl:xh][sel])
        model[xl:xh] = G.predict(x[xl:xh, np.newaxis])
    return model


def subtract_sky(R, sky_sel, args, niter=2, adjustment=None):
    for j in np.arange(niter):
        nwave, nspec = make_avg_spec(R.wave[sky_sel],
                                     safe_division(R.spec[sky_sel],
                                     R.ftf[sky_sel]), binsize=35, knots=None)
        I = interp1d(nwave, nspec, bounds_error=False, kind='quadratic',
                     fill_value='extrapolate')

        R.skysub = R.wave * 0.
        R.sky = R.wave * 0.
        model = R.wave * 0.
        for i in np.arange(R.wave.shape[0]):
            model[i] = I(R.wave[i])
            if adjustment is not None:
                try:
                    sel = np.isfinite(adjustment[i+1])
                    J = interp1d(adjustment[0][sel], adjustment[i+1][sel],
                                 bounds_error=False, kind='quadratic',
                                 fill_value='extrapolate')
                    add = J(R.wave[i])
                except:
                    args.log.warning('Adjustment failed for %s on fiber %i' %
                                     (R.side, i))
                    add = 0.0
                    
            else:
                add = 0.0
            R.sky[i] = model[i] * (R.ftf[i] + add)
            R.skysub[i] = R.spec[i] - R.sky[i]
        residual = safe_division(R.skysub, model)
        sky_sel1 = sky_sel * (np.nanmedian(R.ftf, axis=1) > .5)
        cont = smooth_fiber(residual[:,400:-400], ~sky_sel1, R.wave.shape[0] / 2)
        R.ftf = R.ftf + cont[:, np.newaxis]
        args.log.info('Fiber to Fiber offsets')
        T = Table([R.ifux, R.ifuy, cont], names=['x', 'y', 'offset'])
        args.log.info(T)
    R.skysub = safe_division(R.skysub, R.ftf)
    return R


def extract_source(R, side, lims2, loc, fibinds, args):
    R.skynorm = safe_division(R.sky, R.ftf)
    T = Table.read(op.join(DIRNAME, 'lrs2_config', 'response_%s.dat' % side),
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
                                   lims2, fac=1.0)
    rect_wave, rect_sky = rectify(np.array(R.wave, dtype='float64'),
                                  np.array(R.slam, dtype='float64'),
                                  lims2, fac=1.0)
    R.rect_wave = rect_wave * 1.
    noise = biweight_midvariance(rect_spec, axis=(0,))
    frac, R.flux, R.skyflux, R.fluxerr, R.wei = flux_correction(rect_wave,
                                                         loc, R, fibinds, side,
                                                         rect_spec, rect_sky,
                                                         noise)
    args.log.info('Number of fibers used in extraction: %i' % np.sum(fibinds))
    args.log.info('Fraction of source flux covered in extracted fibers: %0.2f' %
                  np.median(frac))
    return R


def flux_correction(wave, loc, P, inds, side, rect_spec, rect_sky, noise):
    dar_table = Table.read(op.join(DIRNAME, 'lrs2_config', 'dar_%s.dat' % side),
                           format='ascii.fixed_width_two_line')
    X = interp1d(dar_table['wave'], dar_table['x_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    Y = interp1d(dar_table['wave'], dar_table['y_0'], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    frac = wave * 0.
    SF = wave * 0.
    SS = wave * 0.
    N = wave * 0.
    F = rect_spec * 0.
    NX, NY = build_big_fiber_array(P)
    PSF = Moffat2D(amplitude=1., x_0=0., y_0=0., alpha=loc[3], gamma=loc[4])
    for i in np.arange(len(wave)):
        inds1 = np.array(inds*1., dtype=bool)
        x = loc[1] + X(wave[i]) - X(loc[0])
        y = loc[2] + Y(wave[i]) - Y(loc[0])
        PSF.x_0.value = x
        PSF.y_0.value = y
        sind = np.argsort(rect_spec[:, i])
        outlier = ((np.diff(np.sort(rect_spec[:, i])) /
                    np.percentile(rect_spec[:, i], 98)) > 1.)
        for cnt in np.arange(len(outlier)):
            if outlier[cnt]:
                inds1[sind[cnt]] = False
            else:
                break
        for cnt in np.arange(len(outlier))[::-1]:
            if outlier[cnt]:
                inds1[sind[cnt+1]] = False
            else:
                break
        total = PSF(NX, NY).sum()
        wei = PSF(P.ifux[inds1], P.ifuy[inds1])
        frac[i] = wei.sum() / total
        nwei = wei / wei.sum()
        SF[i] = (rect_spec[inds1, i] * nwei).sum() / (nwei**2).sum() / frac[i]
        SS[i] = (rect_sky[inds1, i] * nwei).sum() / (nwei**2).sum() / frac[i]
        N[i] = (np.sqrt((noise[i]**2 * (nwei / (nwei**2).sum())**2).sum()) /
                frac[i])
        F[inds1, i] = nwei
    return frac, SF, SS, N, F


def get_lims(side):
    if side == 'RR':
        lims = [8225, 10565]
        lims2 = [8275., 10500.]
    if side == 'RL':
        lims = [6425, 8460]
        lims2 = [6450., 8400.]
    if side == 'BR':
        lims = [4520, 7010]
        lims2 = [4655., 6960.]
    if side == 'BL':
        lims = [3625, 4670]
        lims2 = [3640., 4640.]
    return lims, lims2


def write_spectrum_out(R):
    names = ['wavelength', 'F_lambda', 'e_F_lambda', 'Sky_lambda']
    hdu = fits.PrimaryHDU(np.array([R.rect_wave, R.flux, R.fluxerr, R.skyflux],
                                   dtype='float32'))
    hdu.header['DWAVE'] = R.rect_wave[1] - R.rect_wave[0]
    hdu.header['WAVE0'] = R.rect_wave[0]
    hdu.header['WAVESOL'] = 'WAVE0 + DWAVE * linspace(0, NAXIS1)'
    hdu.header['WAVEUNIT'] = 'A'
    hdu.header['FLUXUNIT'] = 'ergs/s/cm2/A'
    for i, name in enumerate(names):
        hdu.header['ROW%i' % (i+1)] = name
    for key in R.header.keys():
        if key in hdu.header:
            continue
        hdu.header[key] = R.header[key]
    hdu.writeto(op.join(R.path, 'spectrum_%s.fits' % R.side_name),
                overwrite=True)


def quick_exam(R, nwavebins, lims, side, args, name):
    rect_wave, rect_spec = rectify(np.array(R.wave, dtype='float64'),
                                   np.array(safe_division(R.oldspec,
                                                          R.ftf),
                                            dtype='float64'), lims,
                                   fac=2.5)
    sel = (rect_spec == -999.).sum(axis=0) < 1
    rect_wave = rect_wave[sel]
    rect_spec = rect_spec[:, sel]
    if args.emission:
        rect_spec1 = rect_spec * 1.
        sky = biweight_location(rect_spec, axis=(0, ))
        rect_spec1 = rect_spec - sky[np.newaxis, :]
        rect_spec, mask = convolve_spatially(R.ifux, R.ifuy, rect_spec1, rect_wave,
                                       name, sig_wave=2.5*1.5)
        noise = biweight_midvariance(rect_spec, axis=(0,))
        noise = np.nanmax([0.1 * np.percentile(noise, 95)*np.ones(noise.shape), noise], axis=0)
        ind = np.unravel_index(np.argmax(rect_spec / noise,  axis=None), rect_spec.shape)
        if (args.wave_extract is not None) and (args.extract_side == name):
            cw = np.searchsorted(rect_wave, args.wave_extract)
            sn_image = rect_spec[:, cw]
        else:
            cw = ind[1]
            sn_image = rect_spec[:, cw]
        back = 0.
        wv = rect_wave[cw]
        fits.PrimaryHDU(rect_spec1).writeto('test1.fits', overwrite=True)
        fits.PrimaryHDU(rect_spec / noise).writeto('test2.fits', overwrite=True)
        fits.PrimaryHDU(np.array(mask, dtype=int)).writeto('test3.fits', overwrite=True)

    else:
        func = biweight_location
        Z = np.ma.array(rect_spec, mask=(rect_spec == -999.))
        S = [func(chunk, axis=(1,)) - biweight_location(chunk)
             for chunk in np.array_split(Z, nwavebins, axis=1)]
        B = [biweight_location(chunk)
             for chunk in np.array_split(Z, nwavebins, axis=1)]
        N = [biweight_midvariance(s) for s in S]
        S, B, N = [np.array(v) for v in [S, B, N]]
        SN = [np.percentile(s / n, 99) for s, n in zip(S, N)]
        v = np.argmax(SN)
        sn_image = S[v]
        back = B[v]
        wv = np.median(np.array_split(rect_wave, nwavebins)[v])
    xc, yc, a, g, f, sign, dthresh = find_centroid(sn_image, R.ifux, R.ifuy,
                                                   back)
    xgrid, ygrid, zimage = make_frame(R.ifux, R.ifuy, sn_image)
    mask = mask_sources(xgrid, ygrid, R.ifux, R.ifuy, zimage)
    good = np.setdiff1d(np.arange(rect_spec.shape[0], dtype=int), mask)
    good_mask = np.zeros((rect_spec.shape[0],))
    good_mask[good] = 1.
    good_mask = np.array(good_mask, dtype=bool)
    make_plot(zimage, xgrid, ygrid, R.ifux, R.ifuy, good_mask, R.path,
              side)
    return wv, good_mask, xc, yc, a, g, sign, dthresh


def get_twi_ftf(wave, twi):
    ftf_twi = wave * 0.
    T = np.ma.array(twi, mask=((twi < 1) + np.isnan(twi)))
    fac = biweight_location(T, axis=(1,))[:, np.newaxis] / biweight_location(T)
    ftf_twi = fac * np.ones((twi.shape[1],))
    for i in np.arange(2):
        V =  (twi / ftf_twi).ravel()
        W = wave.ravel()
        v, vind = np.unique(W, return_index=True)
        W = W[vind]
        V = V[vind]
        ind = np.argsort(W)
        sel = np.isfinite(V[ind])
        nsmooth = savgol_filter(V[ind][sel], 35, 1)
        #nwave, smooth = make_avg_spec(wave, twi / ftf_twi)
        I = interp1d(W[ind][sel], nsmooth, kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')
        for i in np.arange(wave.shape[0]):
            ftf_twi[i] = twi[i] / I(wave[i])
            ftf_twi[i] = savgol_filter(ftf_twi[i], 25, 3)
#            N = len(ftf_twi[i])
#            wchunks = np.array_split(wave[i], N / 15)
#            schunks = np.array_split(ftf_twi[i], N / 15)
#            nw = np.array([np.mean(chunk) for chunk in wchunks])
#            sm = np.array([biweight_location(chunk) for chunk in schunks])
#            J = interp1d(nw, sm, kind='quadratic', bounds_error=False,
#                         fill_value='extrapolate')
#            ftf_twi[i] = J(wave[i])
    return ftf_twi


def generate_sky_residual(P, sky_sel, side, lims):
    rect_wave, rect_spec = rectify(np.array(P.wave, dtype='float64'),
                                   np.array(safe_division(P.skysub,
                                                          P.ftf),
                                            dtype='float64'), lims,
                                   fac=1.0)
    nwave, nspec = make_avg_spec(P.wave[sky_sel],
                             safe_division(P.spec[sky_sel],
                             P.ftf[sky_sel]), binsize=35, knots=None)
    btrace = np.zeros((sky_sel.sum(), len(rect_wave)))
    J = interp1d(nwave, nspec, kind='quadratic',
                 fill_value='extrapolate', bounds_error=False)
    S = rect_spec[sky_sel] * 0.
    for kj, ij in enumerate(np.where(sky_sel)[0]):
        I = interp1d(P.wave[ij], P.trace[ij], kind='quadratic',
                     fill_value='extrapolate', bounds_error=False)
        btrace[kj] = I(rect_wave)
        S[kj] = rect_spec[ij] / J(rect_wave)
    s = fits.PrimaryHDU(S)
    t = fits.ImageHDU(btrace)
    v = fits.ImageHDU(rect_wave[:, np.newaxis])
    vv = fits.ImageHDU(np.vstack([nwave, nspec]))
    fits.HDUList([s, t, v, vv]).writeto(op.join(P.path, 'test_%s.fits' % side),
                                        overwrite=True)

def main():
    nwavebins = 20
    min_det_thresh = 5
    seeing = 1.8
    # Load the data
    if args.side == "blue":
        sides = ['BL', 'BR']
        names = ['uv', 'orange']
    else:
        sides = ['RL', 'RR']
        names = ['red', 'farred']
    L = []
    for side, name in zip(sides, names):
        lims, lims2 = get_lims(side)
        R = ReduceLRS2(args.filename, side)

        # Get Mirror Illumination for the given track
        R.get_mirror_illumination()

        # Correct the wavelength solution from sky lines

        # Load the default fiber to fiber and map to each fiber's wavelength
        F = fits.open(op.join(DIRNAME, 'lrs2_config', 'ftf_%s.fits' % side))
        F[0].data = np.array(F[0].data, dtype='float64')
        R.ftf = R.wave * 0.
        for i in np.arange(R.wave.shape[0]):
            I = interp1d(F[0].data[0], F[0].data[i+1], kind='quadratic',
                         bounds_error=False, fill_value=-999.)
            R.ftf[i] = I(R.wave[i])
        R.log.info('Getting fiber to fiber from twilight')
        R.ftf = get_twi_ftf(R.wave, R.twi)
        wv, R.good_mask, xc, yc, a, g, sign, dthresh = quick_exam(R, nwavebins,
                                                                  lims, side,
                                                                  args, name)
        if args.recalculate_wavelength:
            newwave = R.wave * 0.
            args.log.info('Working on the wavelength for side: %s' % side)
            if side in ['BR', 'RL', 'RR']:
                spec = R.oldspec * 1.
#                newwave = get_red_wave(R.wave, R.trace, R.oldspec, R.ftf, R.good_mask,
#                         '%s_skylines.dat' % name, debug=False)
                newwave = fits.open(op.join(DIRNAME, 'lrs2_config', '%s_wavelength.fits' % name))[0].data
                nwave, nspec = make_avg_spec(newwave[R.good_mask],
                                             safe_division(spec, R.ftf)[R.good_mask])
                shift = get_single_shift(nwave, nspec, op.join(DIRNAME, 'lrs2_config', '%s_skylines.dat' % name))
                newwave = newwave + shift
                args.log.info('Shift in wavelength for %s: %0.3f A' % (name, shift))


            else:
                spec = R.twi * 1.
                nwave, nspec = make_avg_spec(R.wave, safe_division(spec, R.ftf))
                newwave = get_new_wave(R.wave, R.trace, spec, R.ftf, R.good_mask,
                                       nwave, nspec)
            wave0 = R.wave * 1.
            R.wave = newwave * 1.
            args.log.info('Max Wave Correction: %0.2f A' %
                          np.max(newwave-wave0))
            args.log.info('Min Wave Correction: %0.2f A' %
                          np.min(newwave-wave0))
            wv, R.good_mask, xc, yc, a, g, sign, dthresh = quick_exam(R, nwavebins, lims, side, args, name)
        L.append([copy(R), sign, xc, yc, wv, dthresh, a, g])
    if (args.wave_extract is not None):
        if (args.extract_side == 'uv') or (args.extract_side == 'red'):
            ind = 0
        if (args.extract_side == 'orange') or (args.extract_side == 'farred'):
            ind = 1
    else:
        ind = np.argmax([l[1] for l in L])
    args.log.info('Point source detected strongest in side: %s' % sides[ind])
    args.log.info('Detection significance: %0.2f' % L[ind][1])
    args.log.info('Detection found at: %0.2f, %0.2f, %0.2f' % (L[ind][2], L[ind][3], L[ind][4]))
    seeing = np.max([1.3, np.min([L[ind][7] * 0.935, 2.5])])
    args.log.info('FWHM: %0.2f' % (seeing))

    if args.wave_extract is not None:
        min_det_thresh = 1.
    if L[ind][1] > min_det_thresh:
        otherind = 1 - ind
        xother, yother = get_x_y_lambda(ind, otherind, L[ind][4],
                                        L[otherind][4], L[ind][2], L[ind][3],
                                        sides)
        L[otherind][2], L[otherind][3] = (xother * 1., yother * 1.)
        for i in np.arange(5, 8):
            L[otherind][i] = L[ind][i]
        for l, side, name in zip(L, sides, names):
            P = l[0]
            d = np.sqrt((l[2] - P.ifux)**2 + (l[3] - P.ifuy)**2)
            sky_sel = d > (seeing * 1.3)
            ext_sel = d < (seeing * 1.3)
            if sky_sel.sum() < 1:
                args.log.warning('Point source found is too bright.')
                args.log.warning('Not enough blank fibers for sky.')
                args.log.warning('Cowardly using 0 for sky.')
                P.sky = P.wave * 0.0
                P.skysub = safe_division(P.spec, P.ftf)
            else:
                if side[0] == 'R':
                    adjustment = np.array(fits.open(op.join(DIRNAME, 'lrs2_config', 'test_%s.fits' % name))[0].data,
                                          dtype='float64')
                else:
                    adjustment = None
                P = subtract_sky(P, sky_sel, args, adjustment=adjustment)
            P.ifupos = np.array([P.ifux, P.ifuy]).swapaxes(0, 1)
            P.skypos = np.array([P.ra, P.dec]).swapaxes(0, 1)

            lims, lims2 = get_lims(side)
            P = extract_source(P, side, lims2, [l[4], l[2], l[3], l[6], l[7]],
                               ext_sel, args)
            P.save(image_list=['image_name', 'error', 'ifupos', 'skypos',
                               'wave', 'twi', 'ftf', 'oldspec', 'sky', 'skysub',
                               'wei', 'trace'],
                   name_list=['image', 'error', 'ifupos', 'skypos', 'wave',
                              'twi', 'ftf', 'spectrum', 'sky', 'skysub', 'weights',
                              'trace'])
            write_spectrum_out(P)
            
            generate_sky_residual(P, sky_sel, side, lims2)
    else:
        args.log.info('No signficant continuum point source found.')
        args.log.info('No spectrum*.fits file will be created.')
        args.log.info('We will assume the field is all sky.')
        args.log.info('A multi*.fits file will still be created.')
        for l, side in zip(L, sides):
            P = l[0]
            lims, lims2 = get_lims(side)
            sky_sel = np.ones(P.ifux.shape, dtype=bool)
            P = subtract_sky(P, sky_sel, args)
            P.ifupos = np.array([P.ifux, P.ifuy]).swapaxes(0, 1)
            P.skypos = np.array([P.ra, P.dec]).swapaxes(0, 1)
            P.save(image_list=['image_name', 'error', 'ifupos', 'skypos',
                               'wave', 'twi', 'ftf', 'oldspec', 'sky', 'skysub',
                               'trace'],
                   name_list=['image', 'error', 'ifupos', 'skypos', 'wave',
                              'twi', 'ftf', 'spectrum', 'sky', 'skysub',
                              'trace'])
            generate_sky_residual(P, sky_sel, side, lims2)

if __name__ == '__main__':
    main()
