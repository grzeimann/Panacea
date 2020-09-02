#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:12:52 2019

@author: gregz
"""

import argparse as ap
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
import sys
import warnings

from scipy.interpolate import interp1d, griddata
from math_utils import biweight
from input_utils import setup_logging
from astropy.convolution import Gaussian1DKernel, convolve, interpolate_replace_nans
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.modeling.models import Polynomial2D, Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
from astropy.stats import sigma_clip, sigma_clipped_stats, mad_std
from fiber_utils_remedy import find_peaks, identify_sky_pixels, get_spectra
from sklearn.decomposition import PCA
from astropy.table import Table

def get_fiber_to_fiber(spec, wave):
    aw = wave.ravel()
    As = spec.ravel()
    inds = np.argsort(aw)
    for j in np.arange(3):
        nw = np.array([np.mean(w) for w in np.array_split(aw[inds], 3500)])
        ns = np.array([biweight(s) for s in np.array_split(As[inds], 3500)])
        I = interp1d(nw, ns, kind='quadratic', bounds_error=False, fill_value='extrapolate')
        ftf = spec * 0.
        for i in np.arange(spec.shape[0]):
            ftf[i] = spec[i] / I(wave[i])
        As = (spec / ftf).ravel()
    return ftf


def get_mastersky(spec, ftf, wave, sel=None):
    if sel is None:
        sel = np.ones((spec.shape[0],), dtype=bool)
    aw = wave[sel].ravel()
    As = ((spec/ ftf)[sel]).ravel()
    inds = np.argsort(aw)
    nw = np.array([np.mean(w) for w in np.array_split(aw[inds], 3500)])
    ns = np.array([biweight(s) for s in np.array_split(As[inds], 3500)])
    I = interp1d(nw, ns, kind='quadratic', bounds_error=False, fill_value='extrapolate')
    sky = spec * 0.
    for i in np.arange(spec.shape[0]):
        sky[i] = I(wave[i])
    return sky, I

def get_rolling_mastersky(spec, ftf, wave, sel=None, size=24):
    if sel is None:
        sel = np.ones((spec.shape[0],), dtype=bool)
    fibnum = np.arange(spec.shape[0])
    sky = spec * 0.
    for i in np.arange(spec.shape[0]):
        selk = (np.abs(fibnum - i) <= size/2.)
        sel2 = sel * selk
        print('Working on Fiber %i with %i fibers' % (i+1, sel2.sum()))
        aw = wave[sel2].ravel()
        As = ((spec[sel2]/ ftf[sel2])).ravel()
        inds = np.argsort(aw)
        nw = np.array([np.mean(w) for w in np.array_split(aw[inds], 3500)])
        ns = np.array([biweight(s) for s in np.array_split(As[inds], 3500)])
        good = np.isfinite(ns)
        I = interp1d(nw[good], ns[good], kind='quadratic', bounds_error=False, fill_value='extrapolate')
        sky[i] = I(wave[i])
    return sky


def correct_amplifier_offsets(y, xp, yp, channel, order=1, kernel=12.):
    xc = xp[np.nanargmax(y)]
    yc = yp[np.nanargmax(y)]
    d = np.sqrt((xp-xc)**2 + (yp-yc)**2)
    k = y * 1.
    k[y==0.] = np.nan
    def split_fit(var, ind=140):
        model = k * 0.
        model[:ind] = convolve(k[:ind], Gaussian1DKernel(kernel), boundary='extend')
        model[ind:] = convolve(k[ind:], Gaussian1DKernel(kernel), boundary='extend')
        return model
    for i in np.arange(5.):
        model = split_fit(k)
        mstd = mad_std(k - model, ignore_nan=True)
        k[(k-model)>2.5*mstd] = np.nan
    model = split_fit(k)
    loc = 2.5
    k[y==0.] = np.nan
    k[d < loc+1.0] = np.nan
    model = split_fit(k)
    bad = np.isnan(k)
    good = ~bad
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    D = fitter(Polynomial2D(1), xp[good], yp[good], (y/model)[good])
    try:
        smodel = D[0](xp, yp)
        mask = D[1]
    except:
        smodel = D[1](xp, yp)
        mask = D[0]
    cor = model * smodel 
    bl, bml = biweight(k[:140]-cor[:140], calc_std=True)
    bl, bmr = biweight(k[140:]-cor[140:], calc_std=True)
    if (bml>0.03) or (bmr>0.03):
        print("Cannot Make Correction")
        z = np.ones(y.shape)
        if channel == 'orange':
            z[:140] = 1.025
            z[140:] = 0.975
        return z, k
    return cor/biweight(cor), k


def get_pca_sky_residuals(data, ncomponents=5):
    pca = PCA(n_components=ncomponents)
    H = pca.fit_transform(data)
    A = np.dot(H, pca.components_)
    return pca, A

def get_residual_map(data, pca):
    res = data * 0.
    for i in np.arange(data.shape[1]):
        sel = np.isfinite(data[:, i])
        coeff = np.dot(data[sel, i], pca.components_.T[sel])
        model = np.dot(coeff, pca.components_)
        res[:, i] = model
    return res


def get_arc_pca(arcskysub, good, mask, components=15):
    X = arcskysub
    X[:, ~mask] = 0.
    X[~good] = 0.
    X = X.swapaxes(0, 1)
    pca, A = get_pca_sky_residuals(X, ncomponents=components)
    return pca

def norm_spec_to_per_A(spec, wave):
    dw = np.diff(wave, axis=1)
    dw = np.hstack([dw[:, 0:1], dw])
    return spec / dw

def rectify(skysub, wave, def_wave):
    skysub_rect = np.zeros((skysub.shape[0], len(def_wave)))
    for i in np.arange(skysub.shape[0]):
        skysub_rect[i] = interp1d(wave[i], skysub[i], kind='linear', bounds_error=False,
                               fill_value=0.0)(def_wave)
    return skysub_rect

def get_apcor(Xc, Yc, d, y):
    x = np.linspace(0, 5.5, 13)
    A = x * 0.
    for i in np.arange(len(x)):
        theta = np.random.rand(20000) * 2. * np.pi
        r = np.random.rand(20000)*0.5 + x[i]
        xr = np.cos(theta) * r
        yr = np.sin(theta) * r
        fr = 0.59 / np.sqrt(3.)
        in_footprint = np.zeros(r.shape, dtype=bool)
        sel = (d > (x[i] - fr)) * (d < (x[i]+0.5+fr))
        for xC, yC in zip(Xc[sel], Yc[sel]):
            in_footprint += np.sqrt((xr-xC)**2 + (yr-yC)**2) < fr
        coverage = (in_footprint > 0).sum() / 20000.
        A[i] = coverage
    c = np.interp(d, x, A, right=0.0)
    apcor = np.nansum(y[d<5.]) / np.nansum(y[d<5.]/c[d<5.])
    return apcor

def find_centroid(pos, y, fibarea, fit_param=None):
    d = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    median, std = biweight(y, calc_std=True)
    sel = d<3.
    y = y - np.nanpercentile(y, 25)
    ind = np.nanargmax(y[sel])
    xc, yc = (pos[sel][ind, 0], pos[sel][ind, 1])
    d = np.sqrt((pos[:, 0] - xc)**2 + (pos[:, 1] - yc)**2)
    median, std = biweight(y[d>3.], calc_std=True)
    a = y[sel][ind]
    G = Gaussian2D(x_mean=xc, y_mean=yc, amplitude=a)
    d = np.sqrt((pos[:, 0] - xc)**2 + (pos[:, 1] - yc)**2)
    sel = (d <= 2.0) * np.isfinite(y)
    fit = LevMarLSQFitter()(G, pos[sel, 0], pos[sel, 1], y[sel])
    new_model= np.sqrt(fit(pos[:, 0], pos[:, 1])*y) 
    new_model[np.isnan(new_model)] = 0.0
    fitquality = False
    if np.nanmax(new_model) > 5 * std:
        fitquality = True
    grid_x, grid_y = np.meshgrid(np.linspace(xc-5., xc+5., 101),
                                 np.linspace(yc-5., yc+5., 101))
    norm = np.sum(fit(grid_x.ravel(), grid_y.ravel())) * 0.1**2
    Xc = pos[:, 0] - xc
    Yc = pos[:, 1] - yc
    apcor = get_apcor(Xc, Yc, d, y)
    return fit.x_mean.value, fit.y_mean.value, fitquality, fit, new_model / norm * fibarea, apcor

def fix_centroid(pos, y, fibarea, fit_param=None):
    median, std = biweight(y, calc_std=True)
    y = y - median
    xc, yc, xs, ys, th = fit_param
    G = Gaussian2D(x_mean=xc, y_mean=yc, x_stddev=xs, y_stddev=ys,
                   theta=th, amplitude=1.)
    G.x_mean.fixed = True
    G.y_mean.fixed = True
    G.x_stddev.fixed = True
    G.y_stddev.fixed = True
    G.theta.fixed = True
    d = np.sqrt((pos[:, 0] - xc)**2 + (pos[:, 1] - yc)**2)
    Xc = pos[:, 0] - xc
    Yc = pos[:, 1] - yc
    sel = (d < 2.0) * np.isfinite(y)
    M = G(pos[sel, 0], pos[sel, 1])
    norm = biweight(y[sel] / M)
    G.amplitude.value = norm
    fit = G
    new_model= np.sqrt(fit(pos[:, 0], pos[:, 1])*y) 
    new_model[np.isnan(new_model)] = 0.0
    fitquality = False
    if np.nanmax(new_model) > 5 * std:
        fitquality = True
    grid_x, grid_y = np.meshgrid(np.linspace(xc-5., xc+5., 101),
                                 np.linspace(yc-5., yc+5., 101))
    norm = np.sum(fit(grid_x.ravel(), grid_y.ravel())) * 0.1**2
    apcor = get_apcor(Xc, Yc, d, y)
    return fit.x_mean.value, fit.y_mean.value, fitquality, fit, new_model / norm * fibarea, apcor


def get_standard(objname, commonwave):
    filename = op.join('/Users/gregz/cure/virus_early/virus_config/'
                       'standards',
                       'm' + objname.lower() + '.dat.txt')
    try:
        wave, standardmag = np.loadtxt(filename, usecols=(0, 1),
                                       unpack=True)
        fnu = 10**(0.4 * (-48.6 - standardmag))
        standard_flam = fnu * 2.99792e18 / wave**2
        standard_wave = wave
        flam = np.interp(commonwave, standard_wave, standard_flam)
        return flam
    except:
        return commonwave * 0.

def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

def get_wave_cor(spec, ftf, wave, mastersky, masterwave):
    Y = spec / ftf
    mask, cont = identify_sky_pixels(mastersky)
    std = mad_std((mastersky-cont)[~mask])
    loc, values = find_peaks((mastersky-cont), thresh=100*std)
    waves = np.interp(loc, np.arange(len(masterwave)), masterwave)
    Waves = np.zeros((280, len(waves))) * np.nan
    Norms = np.zeros((280, len(waves))) * np.nan
    for i in np.arange(Y.shape[0]):
        if np.isfinite(Y[i]).sum() > 200:
            mask, cont = identify_sky_pixels(Y[i])
            std = mad_std((Y[i]-cont)[~mask])
            loc, values = find_peaks((Y[i]-cont), thresh=25*std)
            wav = np.interp(loc, np.arange(Y.shape[1]), wave[i])
            ng = np.isfinite(Y[i]-cont)
            I = interp1d(wave[i][ng], (Y[i]-cont)[ng], kind='quadratic', 
                         bounds_error=False, fill_value=0.0)
            for j in np.arange(len(waves)):
                if len(wav):
                    if np.min(np.abs(wav - waves[j])) < 1.5:
                        Waves[i, j] = wav[np.argmin(np.abs(wav - waves[j]))]
                        Norms[i, j] = I(Waves[i, j])
    return Waves, Norms


def extract_columns(model, chunk, mask=None):
    if model.ndim == 1:
        model = model[: , np.newaxis]
    if mask is None:
        mask = np.isfinite(chunk)
    num1 = np.nansum(chunk * model**2 * mask, axis=0)
    num2 = np.nansum(model**3 * mask, axis=0)
    norm = num1 / num2
    return norm

def get_skyline_mask(sky_rect, mlen=3):
    quick_sky = biweight(sky_rect, axis=0)
    mask, cont = identify_sky_pixels(quick_sky)
    std_sky = mad_std((quick_sky-cont)[~mask])
    loc, values = find_peaks((quick_sky-cont), thresh=15*std_sky)
    loc = np.array(np.round(loc), dtype=int)
    loc = loc[(loc>10) * (loc<(len(quick_sky)-10))]
    Marray = sky_rect * 0.
    for i in np.arange(-mlen, mlen+1):
        Marray[:, loc+i] = np.nan
    return Marray

def get_extraction_model(skysub_rect, sky_rect, def_wave, nchunks=15,
                         func=find_centroid, fit_params=None):
    XC, YC, Nmod, w, XS, YS, TH = ([], [], [], [], [], [], [])
    for chunk, schunk, wi in zip(np.array_split(skysub_rect, nchunks, axis=1),
                                 np.array_split(sky_rect, nchunks, axis=1),
                                 np.array_split(def_wave, nchunks)):
        mod = biweight(chunk, axis=1)
        xc, yc, q, fit, nmod, apcor = func(pos, mod, fibarea, fit_param=fit_params)
        if not too_bright:
            model = nmod 
            model = model / np.nansum(model) * apcor
        else:
            model = mod
            model = model / np.nansum(model) * apcor
        if q:
            Nmod.append(model)
            w.append(np.mean(wi))
            XC.append(xc)
            YC.append(yc)
            XS.append(fit.x_stddev.value)
            YS.append(fit.y_stddev.value)
            TH.append(fit.theta.value)
    return [np.array(xi) for xi in [w, XC, YC, XS, YS, TH]], Nmod, skysub_rect, sky_rect


def get_maxsn_y(skysub, sky, wave, def_wave, pos):
    sky_rect = rectify(sky, wave, def_wave)
    skysub_rect = rectify(skysub, wave, def_wave)
    skyline_mask = get_skyline_mask(sky_rect)
    skysub_rect[np.isnan(skyline_mask)] = np.nan
    G1 = Gaussian1DKernel(1.5)
    smooth = skysub_rect * 0.
    for i in np.arange(skysub_rect.shape[0]):
        smooth[i] = convolve(skysub_rect[i], G1, preserve_nan=True)
    x = pos[:, 0]
    y = pos[:, 1]
    D = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    T = smooth * 0.
    for i in np.arange(smooth.shape[1]):
        T[:, i] = np.nansum(smooth[:, i] * D, axis=1)
    loc1, loc2 = np.unravel_index(np.nanargmax(T), T.shape)
    return T[:, loc2], def_wave[loc2]

def get_source(y, std, spec, pos, fibarea, newftf, wave, sky, check=False):
    # =============================================================================
    # Bright Limit
    # =============================================================================
    too_bright = np.nanmax((y-1.)/std) > 100000.
    
    # =============================================================================
    # Get fit to collapsed spectra
    # =============================================================================
    xc, yc, quality_flag, fit, mod, apcor = find_centroid(pos, y, fibarea)
    d = np.sqrt((xp - xc)**2 + (yp -yc)**2)
    sel = d > (np.max(d) - 2.5)
    dum, std = biweight(y[sel], calc_std=True)
    SN = np.nanmax((y-1.)/std)
    if check:
        return SN
    args.log.info('Maximum S/N fiber: %0.2f' % SN)
    if too_bright:
        sky, I = get_mastersky(spec, newftf, wave, sel=sel)
    return xc, yc, quality_flag, fit, mod, apcor, sky, sel, too_bright

def get_skysub(S, sky, err, d, masksky, channel):
    Sky = biweight(S[d > 5.], axis=0)
    sci = biweight(S[d < 1.5], axis=0)
    masksci = sci > Sky*1.2
    masksky = masksky * (~masksci)
    skysub = S - Sky
    totsky = 0. * skysub
    totsky[:] = Sky
    intermediate = skysub * 1.
    intermediate[:, masksky] = np.nan
    G1 = Gaussian1DKernel(5.5)
    for k in np.arange(S.shape[0]):
        intermediate[k] = interpolate_replace_nans(intermediate[k], G1)
    for k in np.arange(S.shape[1]):
        intermediate[:, k] = interpolate_replace_nans(intermediate[:, k], G1)
    if channel == 'farred':
        mask = masksky * True
        mask[1500:] = False
        mask[:50] = False
        pca = PCA(n_components=55)
    else:
        mask = masksky * True
        mask[1900:] = False
        mask[:50] = False
        pca = PCA(n_components=35)
    y = (skysub - intermediate)[:, mask]
    y[np.isnan(y)] = 0.0
    if mask.sum() > 60:
        pca.fit_transform(y.swapaxes(0, 1))
        res = get_residual_map(skysub - intermediate, pca)
        res[:, ~masksky] = 0.0
        skysub[:] -= res
        totsky[:] += res
    skyfibers = d > 2.
    dummy = skysub * 1.
    ll = np.nanpercentile(dummy, 5, axis=0)
    hh = np.nanpercentile(dummy, 95, axis=0)
    bl, norm = biweight(skysub[:, 200:-200] / err[:, 200:-200], calc_std=True)
    err *= norm
    dummy[dummy > 3.*err] = np.nan
    dummy[(dummy<ll) + (dummy>hh)] = np.nan
    dummy[~skyfibers, :] = np.nan
    dummy[:, masksky] = np.nan
    dummy1 = dummy * 1.
    G = Gaussian2DKernel(7.)
    dummy = convolve(dummy, G, boundary='extend')
    while np.isnan(dummy).sum():
        dummy = interpolate_replace_nans(dummy, G)
    skysub[:] -= dummy
    totsky[:] += dummy
    totsky[:, ~masksky] = 0.0
    return skysub, totsky, dummy1

warnings.filterwarnings("ignore")

DIRNAME = get_script_path()

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--directory",
                    help='''base directory for reductions''',
                    type=str, default="")

parser.add_argument("-c", "--caldirectory",
                    help='''cal directory for reductions''',
                    type=str, default="/work/03946/hetdex/maverick/LRS2/CALS")

parser.add_argument("multiname",
                    help='''e.g., multi_20170126_0000011_exp02_orange.fits''',
                    type=str)

parser.add_argument("-dw", "--delta_wavelength",
                    help='''Delta Wavelength in linear units for output''',
                    default=None, type=float)

parser.add_argument("--fit_params", default=None,
                    help='''e.g., "0.0, 0.0, 1.0, 1.0"''',
                    type=str)

args = None
#args = ['dummy', 'multi_20181116_0000010_exp01_red.fits',
#        '-c', '/Users/gregz/cure/panacea', '-d', '/Users/gregz/cure/panacea']
args = parser.parse_args(args=args)
args.log = setup_logging('lrs2_experiment')

channel_dict = {'uv': 'BL', 'orange': 'BR', 'red': 'RL', 'farred': 'RR'}

# =============================================================================
# Load Data (images and spectra)
# =============================================================================
date = args.multiname.split('_')[1]
channel = args.multiname.split('_')[-1][:-5]
calfile = op.join(args.caldirectory, 'cal_%s_%s.fits' % (date, channel))
m = fits.open(op.join(args.directory, args.multiname))
c = fits.open(calfile)
pos = m['fiber_positions'].data
xp, yp = (pos[:, 0], pos[:, 1])
def_wave = m['extracted_spectrum'].data[0]
trace = c['trace'].data
wave = c['wavelength'].data
fltimage = c['masterFlat'].data
arcimage = c['masterarc'].data
image = m['image'].data
cosmics = m['cosmics'].data
fltspec = get_spectra(fltimage, trace)
arcspec = get_spectra(arcimage, trace)
spec, chi2, errspec = get_spectra(image, trace, array_mod=fltimage)
fltspec, arcspec, spec, errspec = [norm_spec_to_per_A(X, wave)
                              for X in [fltspec, arcspec, spec, errspec]]

fibarea = 1. / 2. * np.sqrt(3.) * 0.59**2


# =============================================================================
# Masking high chi2 values (cosmics and defects that don't flatfield)
# =============================================================================
mask = chi2 > 5
N = mask.shape[0] * mask.shape[1] * 1.
args.log.info('%s: %0.3f masked for chi2' % (args.multiname, (mask.sum() / N)))
spec[mask] = np.nan

if channel == 'uv':
    spec[:, 208:211] = np.nan

# =============================================================================
# Getting Fiber to Fiber
# =============================================================================
ftf = get_fiber_to_fiber(fltspec, wave)
goodfibers = biweight(ftf, axis=1) > 0.5

# =============================================================================
# Reducing Arc
# =============================================================================
arcsky, J = get_mastersky(arcspec, ftf, wave)
yarc = biweight(arcspec / ftf / arcsky, axis=1)
arccor, karc = correct_amplifier_offsets(yarc, xp, yp, channel)
arcftf = ftf * arccor[:, np.newaxis]
arcsky, J = get_mastersky(arcspec, arcftf, wave)
arcskysub = arcspec / arcftf - arcsky
arcskysub_rect = rectify(arcskysub, wave, def_wave)
mask, cont = identify_sky_pixels(J(def_wave))
pca = get_arc_pca(arcskysub_rect, goodfibers, mask, components=95)

# =============================================================================
# Reducing Science
# =============================================================================
sky, I = get_mastersky(spec, ftf, wave)
y = biweight((spec / ftf / sky)[:, 400:600], axis=1)
cor, keep = correct_amplifier_offsets(y, xp, yp, channel)
newftf = ftf * cor[:, np.newaxis]
good = biweight(newftf, axis=1) > 0.5
spec[~good] = np.nan
sel = np.isfinite(keep)
sky, I = get_mastersky(spec, newftf, wave, sel=sel)
iskysub = spec / newftf - sky
ynew, wnew = get_maxsn_y(iskysub, sky, wave, def_wave, pos)
ynew += 1.
yold = biweight(spec / newftf / sky, axis=1)
stdnew = mad_std(ynew[sel])
stdold = mad_std(yold[sel])
wold = def_wave[int(len(def_wave)/2.)]
sel = (np.abs(y - 1.) < 3. * stdold) * good

SN = []
cnt = 0
for y, std in zip([ynew, yold], [stdnew, stdold]):
    sn = get_source(y, std, spec, pos, fibarea, newftf, wave, sky, check=True)
    if cnt == 0:
        args.log.info('Emission SN: %0.2f' % sn)
    if cnt == 1:
        args.log.info('Continuum SN: %0.2f' % sn)
    SN.append(sn)
    cnt += 1
loc = np.argmax(SN)

y = [ynew, yold][loc]
std = [stdnew, stdold][loc]
w = [wnew, wold][loc]
xc, yc, quality_flag, fit, mod, apcor, sky, sel, too_bright = get_source(y, std, spec, pos, fibarea, newftf, wave, sky)



# =============================================================================
# Get Wavelength Adjustment
# =============================================================================
if not too_bright:
    Waves, Norms = get_wave_cor(spec, newftf, wave, I(def_wave), def_wave)
    shift = biweight(Waves - biweight(Waves, axis=0), axis=1)
    shift[np.isnan(shift)] = 0.0
    args.log.info('Wavelength Shift: %0.2f' % biweight(shift))

    wave = wave - shift[:, np.newaxis]
    adj = biweight(Norms / biweight(Norms, axis=0), axis=1)
    #newftf = newftf * adj[:, np.newaxis]
    sky, I = get_mastersky(spec, newftf, wave, sel=sel)


# =============================================================================
# Get PCA residual sky
# =============================================================================
skysub = spec / newftf - sky
skysub_rect = rectify(skysub, wave, def_wave)
spec_rect = rectify(spec / newftf, wave, def_wave)
error_rect = rectify(errspec / newftf, wave, def_wave)
sky_rect = rectify(sky, wave, def_wave)

skysub_rect_orig = skysub_rect * 1.
sky_rect_orig = sky_rect * 1.

# =============================================================================
# Get Extraction Model
# =============================================================================
xs = fit.x_stddev.value
ys = fit.y_stddev.value
th = fit.theta.value
darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % channel_dict[channel])
T = Table.read(darfile, format='ascii.fixed_width_two_line')
xdar = np.interp(w, T['wave'], T['x_0'])
ydar = np.interp(w, T['wave'], T['y_0'])
if args.fit_params is not None:
    fit_p_list = [float(i.replace(' ', '')) for i in args.fit_params.split(',')]
    xc, yc, xs, ys = fit_p_list
xoff = biweight(xc - xdar)
yoff = biweight(yc - ydar)
fit_params = [np.interp(def_wave, T['wave'], T['x_0']+xoff),
              np.interp(def_wave, T['wave'], T['y_0']+yoff),
              xs, ys, th]

N = int(len(def_wave) / 25)
inds = np.arange(int(N/2), len(def_wave), N)
apcor = inds * 0.
W = inds * 0.
for j, i in enumerate(inds):
    W[j] = def_wave[i]
    xc = fit_params[0][i]
    yc = fit_params[1][i]
    d = np.sqrt((pos[:, 0] - xc)**2 + (pos[:, 1] - yc)**2)
    Xc = pos[:, 0] - xc
    Yc = pos[:, 1] - yc
    y = Gaussian2D(x_mean=xc, y_mean=yc, x_stddev=fit_params[2],
                   y_stddev=fit_params[3], theta=fit_params[4])(pos[:, 0], pos[:, 1])
    apcor[j] = get_apcor(Xc, Yc, d, y)
try:
    apcor = np.polyval(np.polyfit(W, apcor, 3), def_wave)
except:
    args.log.warning('Aperture Correction failed due to modeling issue')
    apcor = np.ones(def_wave.shape)
args.log.info('%s: %0.2f %0.2f %0.2f %0.2f' % (args.multiname, fit_params[0][1032], fit_params[1][1032], fit_params[2],
                                                  fit_params[3]))
weight = skysub * 0.
for i in np.arange(skysub.shape[1]):
    xc = fit_params[0][i]
    yc = fit_params[1][i]
    y = Gaussian2D(x_mean=xc, y_mean=yc, x_stddev=fit_params[2],
                   y_stddev=fit_params[3], theta=fit_params[4])(pos[:, 0], pos[:, 1])
    weight[:, i] = y / np.sum(y)
d = np.sqrt((pos[:, 0] - fit_params[0][int(len(def_wave)/2)])**2 +
            (pos[:, 1] - fit_params[1][int(len(def_wave)/2)])**2)

skyline_mask = get_skyline_mask(sky_rect, mlen=5)
skysub_rect, totsky, dummy1 = get_skysub(spec_rect, sky, error_rect, d, np.isnan(skyline_mask.sum(axis=0)), channel)    

spec_rect = extract_columns(weight, skysub_rect_orig)


# =============================================================================
# Get Extraction
# =============================================================================
mask = np.isfinite(skysub_rect)
total_cal = (m['extracted_spectrum'].data[-1] /
              m[0].header['EXPTIME'] /  m[0].header['MILLUM'] /
              m[0].header['THROUGHP'])
spectrum = np.nansum(mask * weight * skysub_rect, axis=0) / np.nansum(mask * weight**2, axis=0)
error = np.sqrt(np.nansum(mask * weight * error_rect**2, axis=0) / np.nansum(mask * weight**2, axis=0))
calibrated = spectrum * total_cal / apcor
mask = np.isfinite(sky_rect)
spectrum_sky = np.nansum(mask * weight * sky_rect, axis=0) / np.nansum(mask * weight**2, axis=0)
calibrated_sky = spectrum_sky * total_cal
spectrum_sum = np.nansum(skysub_rect, axis=0)
calibrated_all = spectrum_sum * total_cal
calibrated_ext = spec_rect * total_cal
calibrated_err = error * total_cal /apcor

fits.PrimaryHDU([def_wave, calibrated, calibrated_sky, calibrated_all, calibrated_ext, spectrum_sum,
                 calibrated_err], header=m[0].header).writeto(
                args.multiname.replace('multi', 'spectrum'), overwrite=True)
fits.PrimaryHDU(np.array(skysub_rect, dtype='float32'), header=m[0].header).writeto(args.multiname.replace('multi', 'skysub'),
                                                         overwrite=True)
fits.PrimaryHDU(np.array(mask, dtype='float32'), header=m[0].header).writeto(args.multiname.replace('multi', 'mask'),
                                                         overwrite=True)
fits.PrimaryHDU(np.array(weight, dtype='float32'), header=m[0].header).writeto(args.multiname.replace('multi', 'weight'),
                                                         overwrite=True)