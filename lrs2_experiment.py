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


def correct_amplifier_offsets(y, xp, yp, order=1, kernel=12.):
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
        return np.ones(y.shape), k
    return cor/biweight(cor), k


def get_pca_sky_residuals(data, ncomponents=5):
    pca = PCA(n_components=ncomponents)
    H = pca.fit_transform(data)
    A = np.dot(H, pca.components_)
    return pca, A

def get_residual_map(data, pca, good):
    res = data * 0.
    for i in np.arange(data.shape[1]):
        sel = good * np.isfinite(data[:, i])
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

def get_apcor(Xc, Yc, d):
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
    mean, median, std = sigma_clipped_stats(y, stdfunc=mad_std)
    y = y - median
    grid_x, grid_y = np.meshgrid(np.linspace(-7., 7., (14*5+1)),
                                 np.linspace(-3.5, 3.5, 7*5+1))
    image = griddata(pos[y>0., :2], y[y>0.], (grid_x, grid_y), method='cubic')
    init_d = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    sel = init_d < 3.0
    xc, yc = (pos[sel, 0][np.nanargmax(y[sel])], pos[sel, 1][np.nanargmax(y[sel])])
    d = np.sqrt((grid_x - xc)**2 + (grid_y - yc)**2)
    sel = (d < 2.) * np.isfinite(image)
    xc = np.sum(image[sel] * grid_x[sel]) / np.sum(image[sel])
    yc = np.sum(image[sel] * grid_y[sel]) / np.sum(image[sel])
    a = y[np.nanargmax(y)]
    G = Gaussian2D(x_mean=xc, y_mean=yc, amplitude=a)
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    d = np.sqrt((pos[:, 0] - xc)**2 + (pos[:, 1] - yc)**2)
    Xc = pos[:, 0] - xc
    Yc = pos[:, 1] - yc
    sel = (d <= 2.0) * np.isfinite(y)
    D = fitter(G, pos[sel, 0], pos[sel, 1], y[sel])
    try:
        fit = D[0]
        dummy = fit(pos[:, 0], pos[:, 1])
        mask = D[1]
    except:
        fit = D[1]
        mask = D[0]
    
    new_model= np.sqrt(fit(pos[:, 0], pos[:, 1])*y) 
    new_model[np.isnan(new_model)] = 0.0
    fitquality = False
    if np.nanmax(new_model) > 5 * std:
        fitquality = True
    grid_x, grid_y = np.meshgrid(np.linspace(xc-5., xc+5., 101),
                                 np.linspace(yc-5., yc+5., 101))
    norm = np.sum(fit(grid_x.ravel(), grid_y.ravel())) * 0.1**2
    apcor = get_apcor(Xc, Yc, d)
    return fit.x_mean.value, fit.y_mean.value, fitquality, fit, new_model / norm * fibarea, apcor

def fix_centroid(pos, y, fibarea, fit_param=None):
    mean, median, std = sigma_clipped_stats(y, stdfunc=mad_std)
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
    apcor = get_apcor(Xc, Yc, d)
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

def get_extraction_model(skysub_rect, sky_rect, def_wave, nchunks=15,
                         niter=1, func=find_centroid, fit_params=None):
    XC, YC, Nmod, w, XS, YS, TH = ([], [], [], [], [], [], [])
    skysub_chunks, sky_chunks, spec_chunks = ([], [], [])
    quick_sky = biweight(sky_rect, axis=0)
    mask, cont = identify_sky_pixels(quick_sky)
    std_sky = mad_std((quick_sky-cont)[~mask])
    loc, values = find_peaks((quick_sky-cont), thresh=15*std_sky)
    loc = np.array(np.round(loc), dtype=int)
    loc = loc[(loc>10) * (loc<(len(quick_sky)-10))]
    Marray = skysub_rect * 0.
    for i in np.arange(-6, 7):
        Marray[:, loc+i] = np.nan
    for chunk, schunk, wi, marray in zip(np.array_split(skysub_rect, nchunks, axis=1),
                                         np.array_split(sky_rect, nchunks, axis=1),
                                         np.array_split(def_wave, nchunks),
                                         np.array_split(Marray, nchunks, axis=1)):
        mod = biweight(chunk, axis=1)
        clean_chunk = chunk * 1.
    
        if fit_params is None:
            fit_param = None
        else:
            fit_param = [np.interp(np.mean(wi), def_wave, fit_params[0]),
                         np.interp(np.mean(wi), def_wave, fit_params[1]),
                         fit_params[2],fit_params[3], fit_params[4]]
        for n in np.arange(1, niter + 1):
            xc, yc, q, fit, nmod, apcor = func(pos, mod, fibarea,
                                               fit_param=fit_param)
            if n == niter:
                print('Iteration %i:' % n, xc, yc, q, '%0.3f' % apcor,
                      fit.x_stddev.value, fit.y_stddev.value, fit.theta.value)
            if not too_bright:
                model = nmod 
                model = model / np.nansum(model) * apcor
            else:
                model = mod
                model = model / np.nansum(model) * apcor
            spectra_chunk = extract_columns(model, chunk)
            model_chunk = model[:, np.newaxis] * spectra_chunk[np.newaxis, :]
            goodpca = np.isfinite(chunk).sum(axis=1) > 0.75 * chunk.shape[1]
            res = get_residual_map(chunk-model_chunk, pca, goodpca)
            blank_image = chunk-model_chunk-res
            bl, bm = biweight(blank_image, axis=0, calc_std=True)
            clean_chunk = chunk - res - bl[np.newaxis, :]
            spectra_chunk = extract_columns(model, clean_chunk)
            mod = biweight(clean_chunk / spectra_chunk[np.newaxis, :], axis=1)
        
        spectra_chunk = extract_columns(model, clean_chunk)
        model_chunk = model[:, np.newaxis] * spectra_chunk[np.newaxis, :]
        spectra_chunk = extract_columns(model, clean_chunk)
        schunk = schunk + res
        skysub_chunks.append(clean_chunk)
        sky_chunks.append(schunk)
        spec_chunks.append(spectra_chunk)
        if q:
            Nmod.append(model)
            w.append(np.mean(wi))
            XC.append(xc)
            YC.append(yc)
            XS.append(fit.x_stddev.value)
            YS.append(fit.y_stddev.value)
            TH.append(fit.theta.value)
    skysub_rect, sky_rect, spec_rect = [np.hstack(x) 
                                        for x in
                                        [skysub_chunks, sky_chunks, spec_chunks]]
    w, xc, yc, xs, ys, th = [np.array(x) for x in [w, XC, YC, XS, YS, TH]]
    Nmod = np.array(Nmod)
    return [w, xc, yc, xs, ys, th], Nmod, skysub_rect, sky_rect, spec_rect


warnings.filterwarnings("ignore")

DIRNAME = get_script_path()

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--directory",
                    help='''base directory for reductions''',
                    type=str, default="")

parser.add_argument("-c", "--caldirectory",
                    help='''cal directory for reductions''',
                    type=str, default="/work/03946/hetdex/maverick/LRS2/CALS")

parser.add_argument("galaxyname",  help='''Name of Galaxy''', type=str)

parser.add_argument("multiname",
                    help='''e.g., multi_20170126_0000011_exp02_orange.fits''',
                    type=str)

parser.add_argument("-dw", "--delta_wavelength",
                    help='''Delta Wavelength in linear units for output''',
                    default=None, type=float)

args = None
#args = ['dummy', 'multi_20181213_0000012_exp01_red.fits',
#        '-c', '/Users/gregz/cure/dummy', '-d', '/Users/gregz/cure/dummy/']
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
spec, chi2 = get_spectra(image, trace, array_mod=fltimage)
fltspec, arcspec, spec = [norm_spec_to_per_A(X, wave)
                          for X in [fltspec, arcspec, spec]]

fibarea = 1. / 2. * np.sqrt(3.) * 0.59**2


# =============================================================================
# Masking high chi2 values (cosmics and defects that don't flatfield)
# =============================================================================
mask = chi2 > 5
spec[mask] = np.nan


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
arccor, karc = correct_amplifier_offsets(yarc, xp, yp)
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
y = biweight(spec / ftf / sky, axis=1)
cor, keep = correct_amplifier_offsets(y, xp, yp)
newftf = ftf * 1. # cor[:, np.newaxis]
good = biweight(newftf, axis=1) > 0.5
spec[~good] = np.nan
sel = np.isfinite(keep)
sky, I = get_mastersky(spec, newftf, wave, sel=sel)
y = biweight(spec / newftf / sky, axis=1)
std = mad_std(y[sel])
sel = (np.abs(y - 1.) < 3. * std) * good

# =============================================================================
# Get fit to collapsed spectra
# =============================================================================
xc, yc, quality_flag, fit, mod, apcor = find_centroid(pos, y, fibarea)
if quality_flag:
    d = np.sqrt((xp - xc)**2 + (yp -yc)**2)
    sel = d > (np.max(d) - 2.5)
    dum, std = biweight(y[sel], calc_std=True)
    args.log.info('Maximum S/N fiber: %0.2f' % (np.nanmax((y-1.)/std)))
    if np.nanmax((y-1.)/std) > 500.:
        sky, I = get_mastersky(spec, newftf, wave, sel=sel)

# =============================================================================
# Bright Limit
# =============================================================================
too_bright = np.nanmax((y-1.)/std) > 500.


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
spec_rect = rectify(spec, wave, def_wave)
sky_rect = rectify(sky, wave, def_wave)

skysub_rect_orig = skysub_rect * 1.
sky_rect_orig = sky_rect * 1.

# =============================================================================
# Get Extraction Model
# =============================================================================
info, Nmod, skysub_rect, sky_rect, spec_rect = get_extraction_model(skysub_rect,
                                                                    sky_rect,
                                                                    def_wave)
w, xc, yc, xs, ys, th = info
darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % channel_dict[channel])
T = Table.read(darfile, format='ascii.fixed_width_two_line')
xdar = np.interp(w, T['wave'], T['x_0'])
ydar = np.interp(w, T['wave'], T['y_0'])
xoff = biweight(xc - xdar)
yoff = biweight(yc - ydar)
XS = biweight(xs)
YS = biweight(ys)
TH = biweight(th)
fit_params = [np.interp(def_wave, T['wave'], T['x_0']+xoff),
              np.interp(def_wave, T['wave'], T['y_0']+yoff),
              XS, YS, TH]
info, Nmod, skysub_rect, sky_rect, spec_rect = get_extraction_model(skysub_rect_orig,
                                                                    sky_rect_orig,
                                                                    def_wave,
                                                             func=fix_centroid,
                                                         fit_params=fit_params)
w, xc, yc, xs, ys, th = info

weight = skysub * 0.
for i in np.arange(skysub.shape[0]):
    fsel = np.isfinite(Nmod[:, i])
    if fsel.sum() > 2:
        osel = (def_wave < w[fsel][0]) + (def_wave > w[fsel][-1])
        weight[i] = interp1d(w[fsel], Nmod[fsel, i], kind='quadratic',
                             fill_value='extrapolate')(def_wave)
        weight[i][osel] = interp1d(w[fsel], Nmod[fsel, i], kind='nearest',
                                 fill_value='extrapolate')(def_wave[osel])

#if not too_bright:
#    spec = extract_columns(weight, skysub_rect)
#    model = spec[np.newaxis, :] * weight
#    res = get_residual_map(skysub_rect_orig-model, pca, good)
#    skysub_rect = skysub_rect_orig - res
#    sky_rect = sky_rect_orig + res

fits.PrimaryHDU(weight, header=m[0].header).writeto(args.multiname.replace('multi', 'weight'),
                                                         overwrite=True)
apcor = np.nansum(weight, axis=0)
weight = weight / apcor[np.newaxis, :]

# =============================================================================
# Get Extraction
# =============================================================================
mask = np.isfinite(skysub_rect)
total_cal = (m['extracted_spectrum'].data[-1] /
              m[0].header['EXPTIME'] /  m[0].header['MILLUM'] /
              m[0].header['THROUGHP'])
spectrum = np.nansum(mask * weight * skysub_rect, axis=0) / np.nansum(mask * weight**2, axis=0)
calibrated = spectrum * total_cal / apcor
mask = np.isfinite(sky_rect)
spectrum_sky = np.nansum(mask * weight * sky_rect, axis=0) / np.nansum(mask * weight**2, axis=0)
calibrated_sky = spectrum_sky * total_cal
spectrum_sum = np.nansum(skysub_rect, axis=0)
calibrated_all = spectrum_sum * total_cal
calibrated_ext = spec_rect * total_cal

fits.PrimaryHDU([def_wave, calibrated, calibrated_sky, calibrated_all, calibrated_ext], header=m[0].header).writeto(
                args.multiname.replace('multi', 'spectrum'), overwrite=True)
fits.PrimaryHDU(skysub_rect_orig, header=m[0].header).writeto(args.multiname.replace('multi', 'skysub'),
                                                         overwrite=True)
