# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:38:06 2018

@author: gregz
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.stats import biweight_midvariance


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
