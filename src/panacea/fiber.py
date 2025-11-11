"""Fiber-level extraction utilities.
"""

import logging
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

from .utils import safe_division

log = logging.getLogger(__name__)


def get_spectra(array_flt, array_trace):
    """Extract per-fiber spectra by sampling along trace centers.

    This method follows the concept of "flat-relative optimal extraction"
    (Zechmeister et al. 2014, A&A, 561, A59), in which the science frame is
    divided by a normalized flat before extraction to mitigate pixel-response
    variations and achieve optimal weighting.

    For each fiber, average the values at the floor/ceil of the trace position
    at each column to approximate the flux centered on the trace.

    Args:
        array_flt: 2D image (e.g., flat-fielded image).
        array_trace: Trace center positions (fibers x columns).

    Returns:
        2D array of shape (fibers, columns) with extracted spectra.
    """
    spectrum = array_trace * 0.0
    x = np.arange(array_flt.shape[1])
    for fiber in np.arange(array_trace.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        try:
            spectrum[fiber] = array_flt[indl, x] / 2.0 + array_flt[indh, x] / 2.0
        except Exception:
            log.warning("Index for getting spectrum out of bounds.")
    return spectrum


def weighted_extraction(image, error, flat, trace, cthresh=8.0):
    """Weighted spectral extraction with cosmic-ray rejection.

    Args:
        image: 2D observed image array.
        error: 2D per-pixel error array.
        flat: 2D flat-field image used for normalization.
        trace: 2D trace positions (fibers x columns).
        cthresh: Threshold for cosmic ray detection.

    Returns:
        Tuple of (spectrum, error_spec, C, Y, Fimage):
        - spectrum: Extracted per-fiber spectra.
        - error_spec: Extracted per-fiber error spectra.
        - C: Boolean cosmic-ray mask.
        - Y: Flat-field normalized image.
        - Fimage: Image labeling fiber assignment per pixel (1..Nfibers).
    """
    E = safe_division(error, flat)
    E[E < 1e-8] = 1e9
    Y = safe_division(image, flat)
    nY = Y * 1.0
    C = np.array(Y * 0.0, dtype=bool)
    # Local import to avoid circular import at module load time
    from .ccd import find_cosmics
    for _ in np.arange(1):
        cosmics = find_cosmics(nY, E, trace, thresh=cthresh, ran=1)
        C = C + cosmics

    x = np.arange(trace.shape[1])
    spectrum = 0.0 * trace
    error_spec = 0.0 * trace
    Fimage = image * 0.0
    for fiber in np.arange(trace.shape[0]):
        T = np.zeros((4, trace.shape[1], 4))
        indl = np.floor(trace[fiber]).astype(int) - 1
        flag = False
        for ss, k in enumerate(np.arange(0, 4)):
            try:
                T[0, :, ss] = Y[indl + k, x]
                T[1, :, ss] = 1.0 / E[indl + k, x] ** 2
                T[2, :, ss] = ~C[indl + k, x]
                T[3, :, ss] = error[indl + k, x] ** 2
                Fimage[indl + k, x] = fiber + 1
            except Exception:
                v = indl + k
                sel = np.where((v >= 0) * (v < Y.shape[0]))[0]
                T[0, sel, ss] = Y[v[sel], x[sel]]
                T[1, sel, ss] = 1.0 / E[v[sel], x[sel]] ** 2
                T[2, sel, ss] = ~C[v[sel], x[sel]]
                T[3, sel, ss] = error[v[sel], x[sel]] ** 2
                Fimage[v[sel], x[sel]] = fiber + 1
                flag = True
        if flag:
            k = 2 if np.mean(indl) > (Y.shape[0] / 2.0) else -1
            v = indl + k
            sel = np.where((v >= 0) * (v < len(x)))[0]
            a = np.sum(T[0, sel] * T[1, sel] * T[2, sel], axis=1)
            b = np.sum(T[1, sel] * T[2, sel], axis=1)
            spectrum[fiber, sel] = safe_division(a, b)
            sel1 = T[2, sel].sum(axis=1) < 2.0
            spectrum[fiber][sel][sel1] = 0.0
            error_spec[fiber][sel][sel1] = 0.0
        else:
            a = np.sum(T[0] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            spectrum[fiber] = safe_division(a, b)
            a = np.sum(T[3] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            error_spec[fiber] = np.sqrt(safe_division(a, b))
            sel = T[2].sum(axis=1) < 2.0
            spectrum[fiber][sel] = 0.0
            error_spec[fiber][sel] = 0.0
    return spectrum, error_spec, C, Y, Fimage


def modify_spectrum(spectrum, error, w):
    """Normalize and fill spectra using wavelength differentials.

    For each fiber, divides spectrum by median dv = diff(w) and fills zeros
    via quadratic interpolation. Error is scaled similarly. Zero-error regions
    are expanded by a few pixels to stabilize extraction.

    Args:
        spectrum: 2D array (fibers x columns) spectral flux.
        error: 2D array matching ``spectrum`` with errors.
        w: 2D array (fibers x columns) of wavelengths.

    Returns:
        Tuple of (spectrum, error) modified in place.
    """
    dw = np.median(np.diff(w, axis=1), axis=0)
    dw = np.hstack([dw[0], dw])
    for i in np.arange(spectrum.shape[0]):
        sel = spectrum[i] == 0.0
        I = interp1d(w[i][~sel], spectrum[i][~sel], kind="quadratic", fill_value="extrapolate")
        spectrum[i] = I(w[i]) / dw
        error[i] /= dw
    bad = error == 0.0
    for i in np.arange(1, 3):
        bad[:, :-i] += bad[:, i:]
        bad[:, i:] += bad[:, :-i]
    error[bad] = 0.0
    return spectrum, error


def correct_ftf(rect, error):
    """Fiber-to-fiber correction for rectified spectra.

    Computes a smooth sky model per column and divides both rect and error by
    the normalized model, unless the object is too bright (then returns inputs).

    Args:
        rect: 2D rectified spectra (fibers x columns).
        error: 2D error array aligned with ``rect``.

    Returns:
        Tuple of (rect, error) corrected, or original arrays if skipped.
    """


    def outlier(y, y1, oi):
        m = np.abs(y[oi] - y1[oi])
        o = (y - y1) > 3.0 * np.median(m)
        return o

    def fit_sky_col(x, y):
        o = y == 0.0
        low = np.percentile(y[~o], 16)
        mid = np.percentile(y[~o], 50)
        high = np.percentile(y[~o], 84)
        flag = False
        if (high - mid) > 2.0 * (mid - low):
            y1 = np.ones(x[~o].shape) * np.percentile(y[~o], 5)
            log.info("Object is too bright for fiber to fiber correction.")
            flag = True
        else:
            y1 = savgol_filter(y[~o], 31, 1)
        I = interp1d(x[~o], y1, kind="quadratic", fill_value="extrapolate")
        y1 = I(x)
        for _ in np.arange(3):
            o += outlier(y, y1, ~o)
            y1 = savgol_filter(y[~o], 51, 1)
            I = interp1d(x[~o], y1, kind="quadratic", fill_value="extrapolate")
            y1 = I(x)
        return y1, o, flag

    x = np.arange(rect.shape[0])
    y = np.median(rect, axis=1)
    f, o, flag = fit_sky_col(x, y)
    ftf = f / np.median(f)
    if not flag:
        return rect / ftf[:, np.newaxis], error / ftf[:, np.newaxis]
    else:
        return rect, error
