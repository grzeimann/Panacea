"""General utilities for Panacea.

Subset migrated from full_lrs2_reduction.py per function_map.md.
"""
from __future__ import annotations

import os.path as op
import tarfile

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.ndimage import percentile_filter

from .wavelength import find_peaks
from .io import get_tarname_from_filename
from .routine import get_mirror_illumination_guider, get_throughput


def power_law(x, c1, c2=0.5, c3=0.15, c4=1.0, sig=2.5):
    """Compute a parametric power-law transform.

    This mirrors the legacy helper used in several modeling routines.

    Args:
        x: Scalar or array input.
        c1: Scale factor.
        c2: Constant added in the denominator.
        c3: Multiplier on the normalized input.
        c4: Exponent on the normalized input.
        sig: Normalization for the input.

    Returns:
        The transformed value(s).
    """
    return c1 / (c2 + c3 * np.power(abs(x / sig), c4))


def safe_division(num, denom, eps=1e-8, fillval=0.0):
    """Safely divide arrays, guarding against tiny or invalid denominators.

    Mirrors the legacy implementation used in spectrum normalization.

    Args:
        num: Numerator array.
        denom: Denominator array (broadcast-compatible with ``num``).
        eps: Values of ``denom`` with absolute value <= eps are considered invalid.
        fillval: Replacement value for invalid divisions.

    Returns:
        Array of the same shape as ``num`` with safe division applied.
    """
    good = np.isfinite(denom) * (np.abs(denom) > eps)
    div = num * 0.0
    if getattr(num, "ndim", 0) == getattr(denom, "ndim", 0):
        div[good] = num[good] / denom[good]
        div[~good] = fillval
    else:
        div[:, good] = num[:, good] / denom[good]
        div[:, ~good] = fillval
    return div



def find_lines(spectrum, trace, nlines, thresh, fib, side=None):
    """Identify arc lines across fibers and align to reference line list.

    This mirrors the legacy heuristic matching between detected peaks in the
    arc spectra and a reference line list table (with columns 'col1', 'col2',
    'col3', 'col4'). It seeds with the brightest expected line and then walks
    through the list by proximity and relative intensity.

    Args:
        spectrum: 2D array (fibers x columns) arc spectra.
        trace: 2D array (fibers x columns) trace positions (for dimensions).
        nlines: astropy.table.Table of reference lines.
        thresh: Peak detection S/N threshold.
        fib: Seed fiber index to start the matching from.
        side: Optional channel name for special-case handling.

    Returns:
        2D array of shape (fibers, n_lines) with detected column positions (0 if missing).
    """
    cont = percentile_filter(spectrum, 15, (1, 101))
    spectrum = spectrum - cont
    loc = []
    ph, pr = ([], [])
    lines = Table(nlines)
    for i, spec in enumerate(spectrum):
        px, ps, py = find_peaks(spec, thresh=thresh)
        sel = np.abs(px - 1032.0) > 0.0
        loc.append(px[sel])
        ph.append(ps[sel])
        pr.append(py[sel])

    if side == 'orange':
        names = ['Hg', 'Cd']
        v = []
        for name in names:
            selhg = lines['col4'] == name
            ma = np.argmax(lines['col3'][selhg])
            sel = np.abs(loc[fib] - lines['col2'][selhg][ma]) < 50.0
            v1 = np.max(pr[fib][sel]) if np.any(sel) else 0.0
            v2 = lines['col3'][selhg][ma]
            v.append([v1, v2])
        if v[0][0] > v[1][0]:
            selhg = lines['col4'] == 'Hg'
            ma = np.argmax(lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            lines['col3'][selhg] *= 1.0 / mxv
            selhg = (lines['col4'] == 'Cd') + (lines['col4'] == 'Ar')
            ma = np.argmax(lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            nrt = v[1][0] / max(v[0][0], 1e-6)
            lines['col3'][selhg] *= nrt / mxv
        else:
            selhg = (lines['col4'] == 'Cd') + (lines['col4'] == 'Ar')
            ma = np.argmax(lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            lines['col3'][selhg] *= 1.0 / mxv
            selhg = lines['col4'] == 'Hg'
            ma = np.argmax(lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            nrt = v[0][0] / max(v[1][0], 1e-6)
            lines['col3'][selhg] *= nrt / mxv

    found_lines = np.zeros((trace.shape[0], len(lines)))
    ls = np.argsort(lines['col3'])[::-1]
    distthresh = 15.0 if side == 'farred' else 50.0
    sel = np.abs(loc[fib] - lines['col2'][ls[0]]) < distthresh
    if np.any(sel):
        ind = np.where(sel)[0][np.argmax(pr[fib][sel])]
        off = loc[fib][ind] - lines['col2'][ls[0]]
        found_lines[fib, ls[0]] = loc[fib][ind]
    else:
        off = 0.0
        found_lines[fib, ls[0]] = 0.0
    y = lines['col2'] + off
    s = np.zeros_like(y)
    pp = np.zeros_like(y)
    pph = np.zeros_like(y)
    s[ls[0]] = 1.0
    pp[ls[0]] = 0.0
    for l in ls[1:]:
        guess = y[l]
        v = np.abs(guess - loc[fib])
        ER = lines['col3'][l] / lines['col3'][ls[0]]
        MR = pr[fib] / max(pr[fib][ind] if 'ind' in locals() else 1.0, 1e-6)
        EE = MR * np.sqrt(1.0 / np.array(ph[fib]) ** 2 + 1.0 / (ph[fib][ind] if 'ind' in locals() else 1.0))
        EE = np.maximum(EE, 0.1 * MR)
        EE = np.maximum(EE, 0.001 * np.ones_like(MR))
        dist = v / 2.0 + np.abs(ER - MR) / EE / 2.0
        if np.min(dist) < 10.0:
            ind1 = np.argmin(dist)
            found_lines[fib, l] = loc[fib][ind1]
            ll = np.where(found_lines[fib] > 0.0)[0][0]
            lh = np.where(found_lines[fib] > 0.0)[0][-1]
            diff0 = [found_lines[fib, ll] - lines['col2'][ll], found_lines[fib, lh] - lines['col2'][lh]]
            m = ((diff0[1] - diff0[0]) / (lines['col2'][lh] - lines['col2'][ll]))
            y = np.array(m * (lines['col2'] - lines['col2'][ll]) + diff0[0] + lines['col2'])
            s[l] = MR[ind1]
            pp[l] = dist[ind1]
            pph[l] = ph[fib][ind1]
    inds = np.where(found_lines[fib] > 0.0)[0]
    delv = []
    for ind2 in inds:
        sel2 = np.where(found_lines[fib, ind2] == found_lines[fib, inds])[0]
        if len(sel2) > 1:
            if np.any(pp[ind2] > pp[inds[sel2]]):
                delv.append(ind2)
                found_lines[fib, ind2] = 0.0
    inds = np.delete(inds, delv)

    for i, line in enumerate(lines):
        if found_lines[fib, i] == 0.0:
            continue
        for j in np.arange(0, fib)[::-1]:
            if len(loc[j]) < 1:
                continue
            k = j + 1
            v = found_lines[k, i]
            while (v == 0.0) and (k < trace.shape[0]):
                k += 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.0:
                found_lines[j, i] = loc[j][np.argmin(m)]
        for j in np.arange(fib + 1, trace.shape[0]):
            if len(loc[j]) < 1:
                continue
            k = j - 1
            v = found_lines[k, i]
            while (v == 0.0) and (k >= 0):
                k -= 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.0:
                found_lines[j, i] = loc[j][np.argmin(m)]
    return found_lines



def create_header_objection(wave, image, func=None):
    """Create an ImageHDU-like object for a 2D spectrum image.

    Args:
        wave: 1D wavelength array.
        image: 2D image to store.
        func: HDU class to instantiate (default astropy.io.fits.ImageHDU).

    Returns:
        ImageHDU with basic WCS-like axis keywords.
    """
    if func is None:
        func = fits.ImageHDU
    hdu = func(np.array(image, dtype='float32'))
    hdu.header['CRVAL1'] = wave[0]
    hdu.header['CRVAL2'] = 1
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CDELT2'] = 1
    hdu.header['CDELT1'] = wave[1] - wave[0]
    return hdu


def build_weight_matrix(x, y, sig=1.5):
    """Build spatial Gaussian weight matrix between fibers.

    Args:
        x: 1D x positions of fibers.
        y: 1D y positions of fibers.
        sig: Gaussian sigma in spatial units.

    Returns:
        2D weight matrix normalized by columns.
    """
    import numpy as _np
    d = _np.sqrt((x - x[:, _np.newaxis])**2 + (y - y[:, _np.newaxis])**2)
    G = _np.exp(-0.5 * (d / sig)**2)
    G = G / G.sum(axis=0)[:, _np.newaxis]
    return G.swapaxes(0, 1)


def mask_skylines_cosmics(wave, rect_spec, name, error):
    """Mask skyline regions and known bad/error pixels.

    Args:
        wave: 1D wavelength array.
        rect_spec: 2D spectra (fibers x wavelength-index).
        name: channel name used to load skyline file.
        error: 2D error array.

    Returns:
        Boolean mask where True indicates masked pixels.
    """
    DIRNAME = op.dirname(op.dirname(op.dirname(__file__)))
    mask1 = rect_spec * 0.0
    if op.exists(op.join(DIRNAME, 'lrs2_config', f'{name}_skylines.dat')):
        T = Table.read(op.join(DIRNAME, 'lrs2_config', f'{name}_skylines.dat'), format='ascii.fixed_width_two_line')
        for w in T['wavelength']:
            mask1[:, np.abs(wave - w) < 6.0] = -1.0
    mask2 = rect_spec * 0.0
    mask2[error == 0.0] = -1.0
    mask2[1:, :] += mask2[:-1, :]
    mask2[:-1, :] += mask2[1:, :]
    if name == 'uv':
        mask1[79:83, 979:987] = -1.0
    mask = (mask1 + mask2) < 0
    return mask


def get_all_cosmics(x, y, ispec, error):
    """Detect spatially inconsistent spikes as cosmics using neighbor median.

    Args:
        x, y: fiber positions.
        ispec: 2D spectra in some window.
        error: 2D error array matching ispec.

    Returns:
        Boolean mask where True indicates likely cosmics.
    """
    D = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    ispec[error == 0.0] = 0.0
    T = ispec * 1.0
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    YY = ispec / T
    YY[np.isnan(YY)] = 0.0
    return YY > 0.2


def convolve_spatially(x, y, spec, wave, name, error, ispec, sig_spatial=0.75, sig_wave=1.5):
    """Spatially and spectrally convolve spectra, masking skylines and cosmics.

    Returns a window around the position of maximum S/N.
    """
    W = build_weight_matrix(x, y, sig=sig_spatial)
    D = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    mask = mask_skylines_cosmics(wave, spec, name, error)
    Z = spec * 1.0
    E = error ** 2
    Z[mask] = np.nan
    E[mask] = np.nan
    ispec[mask] = 0.0
    T = ispec * 1.0
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    YY = ispec / T
    YY[np.isnan(YY)] = 0.0
    Z[YY > 0.2] = np.nan
    E[YY > 0.2] = np.nan
    G = Gaussian1DKernel(sig_wave)
    for i in np.arange(spec.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment='fill', fill_value=0.0)
        E[i, :] = convolve(E[i, :], G, nan_treatment='fill', fill_value=0.0)
    Z_copy = Z * 1.0
    E_copy = np.sqrt(E)
    for i in np.arange(spec.shape[1]):
        Z[:, i] = np.dot(Z[:, i], W)
        E[:, i] = np.dot(E[:, i], W)
    E[:] = np.sqrt(E)
    Y = Z * 0.0
    sel = E > 0.0
    Y[sel] = Z[sel] / E[sel]
    Y[~np.isfinite(Y)] = 0.0
    ind = np.unravel_index(np.nanargmax(Y[:, 50:-50], axis=None), Z[:, 50:-50].shape)
    l1 = ind[1] + 50 - 25
    l2 = ind[1] + 50 + 26
    return ind[1] + 50, Z_copy[:, l1:l2], E_copy[:, l1:l2]


def make_frame(xloc, yloc, data, error, wave, dw, Dx, Dy, wstart=5700., wend=5800., scale=0.4, seeing_fac=1.3):
    """Aggregate a median image over a wavelength window using Gaussian seeing.

    Returns (zgrid, zimage, xgrid, ygrid)
    """
    a, b = data.shape
    x = np.arange(xloc.min()-scale, xloc.max()+1*scale, scale)
    y = np.arange(yloc.min()-scale, yloc.max()+1*scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,) + xgrid.shape)
    area = 3. / 4. * np.sqrt(3.) * 0.59**2
    for k in np.arange(b):
        sel = np.isfinite(data[:, k]) * (error[:, k] != 0.)
        D = np.sqrt((xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid)**2 + (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid)**2)
        W = np.exp(-0.5 / (seeing_fac/2.35)**2 * D**2)
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] * W[sel]).sum(axis=0) / W[sel].sum(axis=0) * (scale**2 / area))
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')
    zimage = np.median(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


def get_objects(basefiles, attrs, full=False):
    """Extract selected FITS header attributes for a set of basefiles.

    If full=True, also appends mirror illumination and throughput (guider based).
    """
    s = []
    for fn in basefiles:
        tarname = get_tarname_from_filename(fn)
        t = tarfile.open(tarname, 'r')
        F = fits.open(t.extractfile('/'.join(fn.split('/')[-4:])))
        s.append([])
        for att in attrs:
            s[-1].append(F[0].header[att])
        if full:
            area = get_mirror_illumination_guider(fn, s[-1][1])
            try:
                throughput = get_throughput(fn, s[-1][1])
            except Exception:
                throughput = 1.0
            s[-1].append(area)
            s[-1].append(throughput)
        t.close()
    return s


def truncate_list(lst):
    """Return 5-element representative subset of a list (first, 3 mids, last)."""
    if len(lst) <= 5:
        return lst
    step = (len(lst) - 1) / 4
    indices = [0, round(step), round(2 * step), round(3 * step), len(lst) - 1]
    return [lst[i] for i in indices]


def get_previous_night(daten):
    """Return YYYYMMDD string for the date one day before ``daten``."""
    from datetime import datetime as _dt, timedelta as _td
    daten_ = _dt(int(daten[:4]), int(daten[4:6]), int(daten[6:])) - _td(days=1)
    return f"{daten_.year:04d}{daten_.month:02d}{daten_.day:02d}"
