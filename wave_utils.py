# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:38:06 2018

@author: gregz
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.stats import biweight_midvariance
from skimage.feature import peak_local_max
from scipy.interpolate import interp1d
from fiber_utils import bspline_x0
from utils import biweight_location
from astropy.table import Table
from astropy.modeling.models import Gaussian1D, Polynomial1D
from astropy.modeling.fitting import LevMarLSQFitter
import matplotlib.pyplot as plt


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


def eval_fit(x1, x2):
    return np.sqrt(biweight_midvariance(x1-x2))


def get_translation(wave, spec, cols, rect_wave, smooth,
                    maxmove=4.):
    I = interp1d(rect_wave, smooth, kind='quadratic',
                 bounds_error=False, fill_value=-999.)
    def f(shifts):
        s = []
        for shift in shifts:
            rect_spec = I(wave[cols] + shift)
            ys = np.ma.array(rect_spec, mask=((rect_spec == 0.) +
                                              (rect_spec == -999.)))
            s.append(eval_fit(spec[cols], ys))
        return shifts[np.argmin(s)], s
    niter = 8
    wid = 2
    center = 0.
    for i in np.arange(niter):
        shifts = np.linspace(-2. * wid + center, 2. * wid + center, 5)
        center, s = f(shifts)
        wid /= 2.

    return center


def get_new_wave(wave, trace, spec, ftf, good_mask, nwave, navg,
                 maxmove=4., nchunks=15, debug=False):
    good_sel = np.where(good_mask)[0]
    wchunks = np.array_split(good_sel, 15)
    wi = np.zeros((len(wchunks), nchunks))
    sh = wi * 0.
    tr = wi * 0.
    xl = wi * 0.
    inds = wi * 0.
    X = np.arange(wave.shape[1])
    for j, wchunk in enumerate(wchunks):
        x = wave[wchunk]
        t = trace[wchunk]
        y = safe_division(spec[wchunk], ftf[wchunk])
        x1 = x.ravel()
        y1 = y.ravel()
        ind = np.argsort(x1)
        chunks = np.array_split(ind, nchunks)
        L = len(wchunk) / 2
        for i, chunk in enumerate(chunks):
            sh[j, i] = get_translation(x1, y1, chunk, nwave, navg)
            wi[j, i] = np.mean(x1[chunk])
            xl[j, i] = np.interp(wi[j, i], x[L], X)
            tr[j, i] = np.interp(xl[j, i], X, np.mean(t, axis=0))
            inds[j, i] = wchunk[L]
    inds = np.array(inds, dtype=int)[:, -1]
    newwave = wave * 1.
    for i, ind in enumerate(inds):
        x = np.interp(wi[i], wave[ind], np.arange(wave.shape[1]))
        p0 = np.polyfit(x / (wave.shape[1] * 1.), wi[i]+sh[i], 3)
        newwave[ind] = np.polyval(p0, np.arange(wave.shape[1]) /
                                  (wave.shape[1] * 1.))
    for i in np.arange(newwave.shape[1]):
        x = trace[inds, i] / (newwave.shape[1] * 1.)
        xfit = trace[:, i] / (newwave.shape[1] * 1.)
        y = newwave[inds, i] * 1.
        newwave[:, i] = np.polyval(np.polyfit(x, y, 3), xfit)
    if debug:
        return newwave, wi, sh, tr, inds, xl
    else:
        return newwave


def get_single_shift(nwave, nspec, skyline_file, debug=False):
    T = Table.read(skyline_file, format='ascii.fixed_width_two_line')
    M = Gaussian1D() + Polynomial1D(0)
    fitter = LevMarLSQFitter()
    match = np.zeros((len(T['wavelength']), 3))
    B, c = bspline_x0(T['wavelength'], nknots=7)
    x = nwave*1.
    y = nspec*1.
    locs1 = peak_local_max(y, min_distance=10, indices=True,
                           threshold_rel=0.001)
    I = interp1d(x, y, kind='quadratic', bounds_error=False,
                 fill_value='extrapolate')
    for j, wi in enumerate(T['wavelength']):
        d = np.abs(x[locs1] - wi)
        ind = np.argmin(d)
        di = d[ind]
        if di < 5.:
            match[j, 0] = wi
            match[j, 2] = locs1[ind]
            nw = np.linspace(wi-3., wi+3., 31)
            M.mean_0 = wi
            M.stddev_0 = 1.8
            M.amplitude_0 = y[locs1][ind]
            M.c0_0 = np.min(I(nw))
            fit = fitter(M, nw, I(nw))
            match[j, 1] = fit.mean_0.value * 1.
    sel = match[:, 1] > 0.
    plt.figure()
    return biweight_location(match[sel, 0] - match[sel, 1])


def get_red_wave(wave, trace, spec, ftf, good_mask, skyline_file,
                 debug=False):
    T = Table.read(skyline_file, format='ascii.fixed_width_two_line')
    M = Gaussian1D() + Polynomial1D(0)
    fitter = LevMarLSQFitter()
    ssel = np.where(good_mask * (np.median(ftf, axis=1) > 0.5))[0]
    match = np.zeros((len(ssel), len(T['wavelength']), 3))
    B, c = bspline_x0(T['wavelength'], nknots=7)
    for k, i in enumerate(ssel):
        x = wave[i]
        y = safe_division(spec[i], ftf[i])
        locs1 = peak_local_max(y, min_distance=10, indices=True,
                               threshold_rel=0.001)
        I = interp1d(x, y, kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')
        for j, wi in enumerate(T['wavelength']):
            d = np.abs(x[locs1] - wi)
            ind = np.argmin(d)
            di = d[ind]
            if di < 5.:
                match[k, j, 0] = wi
                match[k, j, 2] = locs1[ind]
                nw = np.linspace(wi-3., wi+3., 31)
                M.mean_0 = wi
                M.stddev_0 = 1.8
                M.amplitude_0 = y[locs1][ind]
                M.c0_0 = np.min(I(nw))
                fit = fitter(M, nw, I(nw))
                match[k, j, 1] = fit.mean_0.value * 1.
    z = match[:, :, 0] - match[:, :, 1]
    z = np.ma.array(z, mask=(match[:, :, 2] == 0.))
    V = np.ma.median(z, axis=0)
    Wii = np.median(match[:, :, 0], axis=0)
    sel = np.where(Wii > 0.)[0]
    sol = np.linalg.lstsq(c[sel], V[sel])[0]
    Y = np.dot(c, sol)
    JJ = interp1d(T['wavelength'], Y, kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')
    shift = np.zeros((len(ssel)))
    newwave = wave * 1.
    newwave1 = wave * 1.
    for k, i in enumerate(ssel):
        sel = np.where(match[k, :, 2] > 0.)[0]
        yp = match[k, :, 0] - match[k, :, 1]
        shift[k] = biweight_location(JJ(T['wavelength'][sel]) -
                                     yp[sel])

    B, c = bspline_x0(np.arange(280), nknots=11)
    nsel = ssel < 140
    shiftall = np.zeros((wave.shape[0],))
    sol = np.linalg.lstsq(c[ssel[nsel]], shift[nsel])[0]
    Y = np.dot(c, sol)
    shiftall[:140] = Y[:140] * 1.
    nsel = (ssel >= 140) * (ssel < 280)
    sol = np.linalg.lstsq(c[ssel[nsel]], shift[nsel])[0]
    Y = np.dot(c, sol)
    shiftall[140:] = Y[140:] * 1.
    B, c = bspline_x0(np.arange(wave.shape[0]), nknots=7)
    matches = np.ones((match.shape[1], wave.shape[0]))*-99.
    for k in np.arange(match.shape[1]):
        if (match[:, k,  2] > 0.).sum() > (len(ssel) / 2.):
            y = match[:, k, 0] - match[:, k, 1] + shiftall[ssel]
            sol = np.linalg.lstsq(c[ssel], y)[0]
            yp = np.dot(c[ssel], sol)
            nsel = np.abs(y - yp) < (5.*np.sqrt(biweight_midvariance(y-yp)))
            sol = np.linalg.lstsq(c[ssel[nsel]], y[nsel])[0]
            Y = np.dot(c, sol)
            matches[k, :] = Y
            if debug and (k%25 == 0):
                plt.figure()
                plt.scatter(ssel, y)
                plt.plot(np.arange(wave.shape[0]), Y, 'r-')
                plt.ylim([-3., 2.])
                plt.show()
                ans = raw_input('quit?')
                if ans == 'q':
                    sys.exit(1)
    B, c = bspline_x0(T['wavelength'], nknots=7)
    for i in np.arange(wave.shape[0]):
        sel = np.abs(matches[:, i])<3.
        sol = np.linalg.lstsq(c[sel], matches[sel, i])[0]
        Y = np.dot(c, sol)
        J = interp1d(T['wavelength'], Y, kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')
        nw = J(wave[i]) + wave[i]
        newwave[i] = wave[i] + J(nw) - shiftall[i]

        if debug and (i%20 == 0):
            plt.figure()
            plt.scatter(T['wavelength'][sel], matches[sel, i])
            plt.plot(wave[i], J(nw))
            plt.ylim([-3., 2.])
            plt.show()
            ans = raw_input('quit?')
            if ans == 'q':
                sys.exit(1)
    for i in np.arange(newwave.shape[1]):
        x = trace[ssel, i] / (newwave.shape[1] * 1.)
        xfit = trace[:, i] / (newwave.shape[1] * 1.)
        y = newwave[ssel, i] * 1.
        newwave1[:, i] = np.polyval(np.polyfit(x, y, 7), xfit)

    return newwave1
