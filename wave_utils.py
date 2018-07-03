# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:38:06 2018

@author: gregz
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.stats import biweight_midvariance


def eval_fit(x1, x2):
    return np.sqrt(biweight_midvariance(x1-x2))


def get_translation(wave, spec, cols, avg, xp, smooth, rect_wave,
                    maxmove=4.):
    dw = np.diff(wave)
    dw = np.hstack([dw[0], dw])
    I = interp1d(wave, spec / dw, kind='quadratic',
                 bounds_error=False, fill_value=-999.)
    def f(shifts):
        s = []
        for shift in shifts:
            rect_spec = I(rect_wave-shift)
            ys = np.ma.array(rect_spec, mask=((rect_spec == 0.) +
                                              (rect_spec == -999.)))
            py1 = ys[~avg.mask][cols, np.newaxis]
            py2 = smooth[~avg.mask][cols, np.newaxis]
            ratio = py1 / py2
            sel = np.squeeze(~ratio.mask)
            norm1 = np.polyval(np.polyfit(xp[sel], ratio[sel], 2), xp)[:,
                                                                    np.newaxis]
            s.append(eval_fit(py1 / norm1, py2))
        return shifts[np.argmin(s)], s
    niter = 8
    wid = 2
    center = 0.
    for i in np.arange(niter):
        shifts = np.linspace(-2. * wid + center, 2. * wid + center, 5)
        center, s = f(shifts)
        wid /= 2.

    return center


def get_new_wave(wave, trace, spec, rect_wave, avg, smooth, maxmove=4.,
                 nchunks=15):
    col_chunk = np.array_split(np.arange(20, (~avg.mask).sum()-20), nchunks)
    x = np.arange((~avg.mask).sum())
    newwave = wave * 1.
    NFIBS = wave.shape[0]
    shifts = []
    warray, xarray, yarray, oarray = ([], [], [], [])
    inds = np.array(np.hstack([np.arange(1, NFIBS, 8), NFIBS-1]), dtype=int)
    for ind in inds:
        for j, cols in enumerate(col_chunk):
            newwave[ind] = wave[ind] * 1.
            xp = x[cols]
            x0 = np.searchsorted(wave[ind],
                                 np.median(rect_wave[cols]))
            yarray.append(trace[ind, x0])
            xarray.append(x0)
            totshift = 0.
            totshift = get_translation(wave[ind], spec[ind], cols, avg, xp,
                                       smooth, rect_wave,
                                       maxmove=maxmove)
            shifts.append(totshift)
            warray.append(wave[ind, x0] + totshift)
            oarray.append(wave[ind, x0])
    x, y, w = [np.hstack(i) for i in [xarray, yarray, warray]]
    wi, sh = [np.array(v).reshape(len(inds), nchunks)
              for v in [oarray, shifts]]
    newwave = wave * 1.
    for i, ind in enumerate(inds):
        x = np.interp(wi[ind], wave[ind], np.arange(wave.shape[1]))
        p0 = np.polyfit(x / (wave.shape[1] * 1.), wi[i]+sh[i], 3)
        newwave[ind] = np.polyval(p0, np.arange(wave.shape[1]) /
                                  (wave.shape[1] * 1.))
    for i in np.arange(newwave.shape[1]):
        x = trace[inds, i] / (newwave.shape[1] * 1.)
        xfit = trace[:, i] / (newwave.shape[1] * 1.)
        y = newwave[inds, i] * 1.
        newwave[:, i] = np.polyval(np.polyfit(x, y, 3), xfit)
    return newwave
