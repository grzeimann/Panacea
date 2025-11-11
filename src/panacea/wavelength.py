"""Wavelength calibration utilities.

Functions migrated per function_map.md.
"""
from __future__ import annotations

import numpy as np
from astropy.stats import biweight_midvariance
import logging
from astropy.convolution import Gaussian1DKernel, convolve
from .utils import find_lines, robust_polyfit

def find_peaks(y, thresh = 8.0):
    """Locate significant peaks in a 1D array using slope changes and S/N.

    Mirrors legacy implementation: detects zero-crossings in the derivative,
    filters by threshold in units of robust std (biweight), and refines peak
    positions using a quadratic interpolation about the discrete maximum.

    Args:
        y: 1D array to search for peaks.
        thresh: Threshold in units of robust sigma for accepting peaks.

    Returns:
        Tuple of (peak_loc, normalized_peaks, peaks):
        - peak_loc: Sub-pixel peak locations.
        - normalized_peaks: Peak heights divided by robust sigma.
        - peaks: Raw peak heights at rounded peak_loc indices.
    """
    def _get_peaks(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.0
        inds[1] = XN + 0.0
        inds[2] = XN + 1.0
        inds = np.array(inds, dtype=int)
        Peaks = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2.0 * (flat[inds[2]] - 2.0 * flat[inds[1]] + flat[inds[0]])))
        return Peaks

    diff_array = y[1:] - y[:-1]
    loc = np.where((diff_array[:-1] > 0.0) * (diff_array[1:] < 0.0))[0]
    peaks = y[loc + 1]
    std = np.sqrt(biweight_midvariance(y))
    loc = loc[peaks > (thresh * std)] + 1
    peak_loc = _get_peaks(y, loc)
    peaks = y[np.round(peak_loc).astype(int)]
    return peak_loc, peaks / std, peaks



def get_wavelength_from_arc(image, trace, lines, side, amp, date, otherimage=None):
    """Derive per-fiber wavelength solutions from arc images.

    Mirrors legacy logic: extract per-fiber spectra from the arc image, detect
    emission lines, map observed pixel positions to reference wavelengths, and
    fit a 3rd-order polynomial dispersion per fiber. For the 'uv' channel after
    2016-11-01, uses a smoothed median across neighboring fibers and optionally
    a reference arc image to stabilize line detection.

    Args:
        image: 2D arc image (pixels_y x pixels_x).
        trace: 2D array (fibers x columns) of trace center positions.
        lines: astropy.table.Table with columns 'col1' (wavelength) and
            'col2' (expected column), 'col3' (relative intensity), 'col4' (name).
        side: Channel name string (e.g., 'uv', 'orange', 'red', 'farred').
        amp: Amplifier identifier (unused; for logging/compatibility).
        date: Integer-like yyyymmdd for special-case behavior.
        otherimage: Optional second arc image for 'uv' stabilization.

    Returns:
        2D array same shape as ``trace`` containing wavelength for each pixel
        along columns per fiber.
    """


    log = logging.getLogger(__name__)

    # Local import to avoid circular import: wavelength -> fiber -> ccd -> sky -> wavelength
    from .fiber import get_spectra
    spectrum = get_spectra(image, trace)
    fib = int(np.argmax(np.median(spectrum, axis=1)))
    if side == 'uv' and int(date) > 20161101:
        thresh = 5.0
        spectrum2 = spectrum * 0.0
        for i in np.arange(trace.shape[0]):
            ll = int(np.max([0, i - 4]))
            hl = int(np.min([trace.shape[0], i + 5]))
            spectrum2[i] = np.median(spectrum[ll:hl], axis=0)
            G = Gaussian1DKernel(1.5)
            spectrum2[i] = convolve(spectrum2[i], G)
        spectrum = spectrum2 * 1.0
        if otherimage is not None:
            spectrum1 = get_spectra(otherimage, trace)
            spectrum2 = spectrum1 * 0.0
            for i in np.arange(trace.shape[0]):
                ll = int(np.max([0, i - 4]))
                hl = int(np.min([trace.shape[0], i + 5]))
                spectrum2[i] = np.median(spectrum1[ll:hl], axis=0)
                G = Gaussian1DKernel(1.5)
                spectrum2[i] = convolve(spectrum2[i], G)
            spectrum1 = spectrum2 * 1.0
            fib = int(trace.shape[0] / 2)
            found_lines = find_lines(spectrum, trace, lines, thresh, fib)
            found_lines1 = find_lines(spectrum1, trace, lines, thresh, fib)
            sel = (found_lines > 0.0) * (found_lines1 > 0.0)
            for i in np.arange(found_lines.shape[0]):
                found_lines[i] = (
                    found_lines1[i]
                    + np.median(found_lines[i][sel[i]] - found_lines1[i][sel[i]])
                )
                found_lines[i][found_lines1[i] == 0.0] = 0.0
        else:
            found_lines = find_lines(spectrum, trace, lines, thresh, fib)
    else:
        thresh = 3.0
        found_lines = find_lines(spectrum, trace, lines, thresh, fib, side)

    x = np.arange(trace.shape[1])
    found_lines[found_lines > 2060] = 0.0
    found_lines1 = found_lines * 1.0

    qft = np.zeros((len(lines),))
    for i, line in enumerate(lines):
        if np.sum(found_lines[:, i]) < (0.5 * trace.shape[0]):
            found_lines[:, i] = 0.0
            continue
        ind = np.array(found_lines[:, i], dtype=int)
        xt = trace[np.arange(trace.shape[0]), ind]
        yt = robust_polyfit(xt, found_lines[:, i])
        sel = found_lines1[:, i] > 0.0
        qft[i] = np.std(found_lines1[sel, i] - yt[sel])
        if qft[i] > 0.25:
            found_lines[:, i] = 0.0
            continue
        found_lines[:, i] = yt

    wave = trace * 0.0
    res = np.zeros((trace.shape[0],))
    for j in np.arange(trace.shape[0]):
        sel = found_lines[j, :] > 0.0
        if sel.sum() < 3:
            continue
        wave[j] = np.polyval(np.polyfit(found_lines[j, sel], lines['col1'][sel], 3), x)
        res[j] = np.std(np.interp(found_lines[j, sel], x, wave[j]) - lines['col1'][sel])

    missing = np.where(np.all(wave == 0.0, axis=1))[0]
    if len(missing):
        good = np.where(~np.all(wave == 0.0, axis=1))[0]
        for j in np.arange(trace.shape[1]):
            xx = trace[good, j]
            yy = wave[good, j]
            wave[missing, j] = np.polyval(np.polyfit(xx, yy, 3), trace[missing, j])

    log.info('Min, Max Wave: %0.2f, %0.2f', wave.min(), wave.max())
    log.info('Mean Res, Median Res: %0.3f, %0.3f', np.mean(res), np.median(res))
    return wave
