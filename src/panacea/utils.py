"""General utilities for Panacea.

Subset migrated from full_lrs2_reduction.py per function_map.md.
"""

from __future__ import annotations

import os.path as op
import tarfile
import sys
from datetime import datetime, timedelta
from importlib import resources

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.ndimage import percentile_filter

from .io import get_tarname_from_filename

# Standard star canonical names used to identify calibration frames
STANDARD_NAMES = [
    "HD_19445",
    "SA95-42",
    "GD50",
    "G191B2B",
    "HILTNER_600",
    "G193-74",
    "PG0823+546",
    "HD_84937",
    "GD108",
    "FEIGE_34",
    "HD93521",
    "GD140",
    "HZ_21",
    "FEIGE_66",
    "FEIGE_67",
    "G60-54",
    "HZ_44",
    "GRW+70_5824",
    "BD+26+2606",
    "BD+33_2642",
    "G138-31",
    "WOLF_1346",
    "BD_+17_4708",
    "FEIGE_110",
    "GD248",
    "HZ_4",
    "BD+40_4032",
    "HILTNER_102",
    "BD_+26_2606",
    "GD_248",
    "FEIGE_56",
    "FEIGE_92",
    "HZ_15",
    "FEIGE_98",
    "BD+08_2015",
    "BD+25_3941",
    "FEIGE_15",
    "FEIGE_25",
    "SA_95-42",
    "BD+28_4211",
    "HR6203",
]


def check_if_standard(objname: str) -> bool:
    """Return True if objname matches a known spectrophotometric standard.

    Case-insensitive substring match against STANDARD_NAMES.
    """
    for standard in STANDARD_NAMES:
        if standard.lower() in objname.lower():
            return True
    return False


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
    from .wavelength import find_peaks

    for i, spec in enumerate(spectrum):
        px, ps, py = find_peaks(spec, thresh=thresh)
        sel = np.abs(px - 1032.0) > 0.0
        loc.append(px[sel])
        ph.append(ps[sel])
        pr.append(py[sel])

    if side == "orange":
        names = ["Hg", "Cd"]
        v = []
        for name in names:
            selhg = lines["col4"] == name
            ma = np.argmax(lines["col3"][selhg])
            sel = np.abs(loc[fib] - lines["col2"][selhg][ma]) < 50.0
            v1 = np.max(pr[fib][sel]) if np.any(sel) else 0.0
            v2 = lines["col3"][selhg][ma]
            v.append([v1, v2])
        if v[0][0] > v[1][0]:
            selhg = lines["col4"] == "Hg"
            ma = np.argmax(lines["col3"][selhg])
            mxv = lines["col3"][selhg][ma]
            lines["col3"][selhg] *= 1.0 / mxv
            selhg = (lines["col4"] == "Cd") + (lines["col4"] == "Ar")
            ma = np.argmax(lines["col3"][selhg])
            mxv = lines["col3"][selhg][ma]
            nrt = v[1][0] / max(v[0][0], 1e-6)
            lines["col3"][selhg] *= nrt / mxv
        else:
            selhg = (lines["col4"] == "Cd") + (lines["col4"] == "Ar")
            ma = np.argmax(lines["col3"][selhg])
            mxv = lines["col3"][selhg][ma]
            lines["col3"][selhg] *= 1.0 / mxv
            selhg = lines["col4"] == "Hg"
            ma = np.argmax(lines["col3"][selhg])
            mxv = lines["col3"][selhg][ma]
            nrt = v[0][0] / max(v[1][0], 1e-6)
            lines["col3"][selhg] *= nrt / mxv

    found_lines = np.zeros((trace.shape[0], len(lines)))
    ls = np.argsort(lines["col3"])[::-1]
    distthresh = 15.0 if side == "farred" else 50.0
    sel = np.abs(loc[fib] - lines["col2"][ls[0]]) < distthresh
    if np.any(sel):
        ind = np.where(sel)[0][np.argmax(pr[fib][sel])]
        off = loc[fib][ind] - lines["col2"][ls[0]]
        found_lines[fib, ls[0]] = loc[fib][ind]
    else:
        off = 0.0
        found_lines[fib, ls[0]] = 0.0
    y = lines["col2"] + off
    s = np.zeros_like(y)
    pp = np.zeros_like(y)
    pph = np.zeros_like(y)
    s[ls[0]] = 1.0
    pp[ls[0]] = 0.0
    for line_idx in ls[1:]:
        guess = y[line_idx]
        v = np.abs(guess - loc[fib])
        ER = lines["col3"][line_idx] / lines["col3"][ls[0]]
        MR = pr[fib] / max(pr[fib][ind] if "ind" in locals() else 1.0, 1e-6)
        EE = MR * np.sqrt(
            1.0 / np.array(ph[fib]) ** 2
            + 1.0 / (ph[fib][ind] if "ind" in locals() else 1.0)
        )
        EE = np.maximum(EE, 0.1 * MR)
        EE = np.maximum(EE, 0.001 * np.ones_like(MR))
        dist = v / 2.0 + np.abs(ER - MR) / EE / 2.0
        if np.min(dist) < 10.0:
            ind1 = np.argmin(dist)
            found_lines[fib, line_idx] = loc[fib][ind1]
            ll = np.where(found_lines[fib] > 0.0)[0][0]
            lh = np.where(found_lines[fib] > 0.0)[0][-1]
            diff0 = [
                found_lines[fib, ll] - lines["col2"][ll],
                found_lines[fib, lh] - lines["col2"][lh],
            ]
            m = (diff0[1] - diff0[0]) / (lines["col2"][lh] - lines["col2"][ll])
            y = np.array(
                m * (lines["col2"] - lines["col2"][ll]) + diff0[0] + lines["col2"]
            )
            s[line_idx] = MR[ind1]
            pp[line_idx] = dist[ind1]
            pph[line_idx] = ph[fib][ind1]
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
    hdu = func(np.array(image, dtype="float32"))
    hdu.header["CRVAL1"] = wave[0]
    hdu.header["CRVAL2"] = 1
    hdu.header["CRPIX1"] = 1
    hdu.header["CRPIX2"] = 1
    hdu.header["CTYPE1"] = "pixel"
    hdu.header["CTYPE2"] = "pixel"
    hdu.header["CDELT2"] = 1
    hdu.header["CDELT1"] = wave[1] - wave[0]
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
    d = np.sqrt((x - x[:, np.newaxis]) ** 2 + (y - y[:, np.newaxis]) ** 2)
    G = np.exp(-0.5 * (d / sig) ** 2)
    G = G / G.sum(axis=0)[:, np.newaxis]
    return G.swapaxes(0, 1)


def get_config_file(filename):
    """Return a Traversable to a config file inside panacea/lrs2_config.

    Uses importlib.resources to locate packaged data regardless of install mode.

    Args:
        filename: Name of the file to locate.

    Returns:
        pathlib.Traversable: Traversable to the file.
    """
    return resources.files("panacea") / "lrs2_config" / filename


def read_arc_lines(file_obj):
    """Read arc line list files robustly from whitespace-separated data.

    The packaged line-list files under panacea/lrs2_config/lines_*.dat are
    whitespace/tab-separated with comment lines starting with '#'. This parser
    reads the first four columns per data line and returns an Astropy Table
    with the expected column names ['col1','col2','col3','col4'].

    Args:
        file_obj: Text file-like object opened for reading.

    Returns:
        astropy.table.Table: Table with columns ['col1','col2','col3','col4']
            corresponding to wavelength, approx_x, relative intensity, and name.
    """
    rows = []
    for raw in file_obj:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        # Require at least 4 columns: wavelength, approx_x, rel_intensity, name
        if len(parts) < 4:
            continue
        try:
            lam = float(parts[0])
            approx_x = float(parts[1])
            rel = float(parts[2])
            name = parts[3]
            rows.append((lam, approx_x, rel, name))
        except Exception:
            # Skip malformed lines gracefully
            continue
    t = Table(rows=rows, names=["col1", "col2", "col3", "col4"])
    # Provide a numpy-like shape attribute for convenience in tests and callers
    # Astropy Table does not define .shape by default, but it allows setting arbitrary attributes.
    try:
        t.shape = (len(t), len(t.colnames))
    except Exception:
        pass
    return t


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
    mask1 = rect_spec * 0.0
    cfg = get_config_file(f"{name}_skylines.dat")
    try:
        if hasattr(cfg, "is_file") and cfg.is_file():
            wavelengths = []
            with cfg.open("r") as f:
                for raw in f:
                    line = raw.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    try:
                        wavelengths.append(float(parts[0]))
                    except Exception:
                        continue
            for w in wavelengths:
                mask1[:, np.abs(wave - w) < 6.0] = -1.0
    except FileNotFoundError:
        pass
    mask2 = rect_spec * 0.0
    mask2[error == 0.0] = -1.0
    mask2[1:, :] += mask2[:-1, :]
    mask2[:-1, :] += mask2[1:, :]
    if name == "uv":
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
    D = np.sqrt((x - x[:, np.newaxis]) ** 2 + (y - y[:, np.newaxis]) ** 2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    ispec[error == 0.0] = 0.0
    T = ispec * 1.0
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    # Guard against division by zero/invalid to avoid RuntimeWarning and spurious inf/NaN
    YY = np.divide(
        ispec, T, out=np.zeros_like(ispec), where=np.isfinite(T) & (np.abs(T) > 0)
    )
    return YY > 0.2


def convolve_spatially(
    x, y, spec, wave, name, error, ispec, sig_spatial=0.75, sig_wave=1.5
):
    """Spatially and spectrally convolve spectra, masking skylines and cosmics.

    Applies a 1D Gaussian convolution along wavelength and a spatial
    Gaussian-weighted sum across neighboring fibers after masking strong
    skylines and likely cosmic rays. Identifies the wavelength index with
    the highest median S/N and returns a small window around it along with
    the corresponding per-fiber convolved signal and error slices.

    Args:
        x (np.ndarray): 1D x positions of fibers (length Nfibers).
        y (np.ndarray): 1D y positions of fibers (length Nfibers).
        spec (np.ndarray): 2D rectified spectra array with shape (Nfibers, Nwave).
        wave (np.ndarray): 1D wavelength array of length Nwave.
        name (str): Spectral channel name used for skyline masking (e.g., 'uv').
        error (np.ndarray): 2D error array aligned with ``spec``.
        ispec (np.ndarray): 2D intermediate spectrum used to detect cosmics
            (same shape as ``spec``).
        sig_spatial (float, optional): Sigma of the spatial Gaussian used to
            build the fiber-to-fiber weights. Defaults to 0.75.
        sig_wave (float, optional): Sigma of the 1D Gaussian kernel applied along
            wavelength. Defaults to 1.5.

    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            - loc: Integer wavelength index of the S/N peak.
            - sdimage: 2D array (Nfibers, 51) of convolved signal centered on
              ``loc``.
            - sderror: 2D array (Nfibers, 51) of corresponding errors.
    """
    W = build_weight_matrix(x, y, sig=sig_spatial)
    D = np.sqrt((x - x[:, np.newaxis]) ** 2 + (y - y[:, np.newaxis]) ** 2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    mask = mask_skylines_cosmics(wave, spec, name, error)
    Z = spec * 1.0
    E = error**2
    Z[mask] = np.nan
    E[mask] = np.nan
    ispec[mask] = 0.0
    T = ispec * 1.0
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    # Guard against division by zero/invalid to avoid RuntimeWarning and spurious inf/NaN
    YY = np.divide(
        ispec, T, out=np.zeros_like(ispec), where=np.isfinite(T) & (np.abs(T) > 0)
    )
    Z[YY > 0.2] = np.nan
    E[YY > 0.2] = np.nan
    G = Gaussian1DKernel(sig_wave)
    for i in np.arange(spec.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment="fill", fill_value=0.0)
        E[i, :] = convolve(E[i, :], G, nan_treatment="fill", fill_value=0.0)
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


def make_frame(
    xloc,
    yloc,
    data,
    error,
    wave,
    dw,
    Dx,
    Dy,
    wstart=5700.0,
    wend=5800.0,
    scale=0.4,
    seeing_fac=1.3,
):
    """Create a collapsed spatial image over a wavelength window.

    Builds per-wavelength model images by distributing fiber fluxes onto a
    regular spatial grid with a circular Gaussian PSF and then takes the
    median across a wavelength interval [wstart, wend].

    Args:
        xloc (np.ndarray): 1D x positions of fibers.
        yloc (np.ndarray): 1D y positions of fibers.
        data (np.ndarray): 2D rectified flux array (Nfibers, Nwave).
        error (np.ndarray): 2D error array aligned with ``data``.
        wave (np.ndarray): 1D wavelength grid of length Nwave.
        dw (np.ndarray): 1D delta-wavelength per column (unused legacy arg).
        Dx (np.ndarray): 1D array of DAR x offsets per wavelength index.
        Dy (np.ndarray): 1D array of DAR y offsets per wavelength index.
        wstart (float, optional): Start wavelength for collapsing. Defaults to 5700.
        wend (float, optional): End wavelength for collapsing. Defaults to 5800.
        scale (float, optional): Spatial grid step size in IFU units. Defaults to 0.4.
        seeing_fac (float, optional): FWHM multiplier controlling the Gaussian PSF
            sigma as seeing_fac/2.35. Defaults to 1.3.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            - zgrid: 3D cube with per-wavelength model images (Nwave, Ny, Nx).
            - zimage: 2D collapsed median image over [wstart, wend].
            - xgrid: 2D x coordinate grid.
            - ygrid: 2D y coordinate grid.
    """
    a, b = data.shape
    x = np.arange(xloc.min() - scale, xloc.max() + 1 * scale, scale)
    y = np.arange(yloc.min() - scale, yloc.max() + 1 * scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,) + xgrid.shape)
    area = 3.0 / 4.0 * np.sqrt(3.0) * 0.59**2
    for k in np.arange(b):
        sel = np.isfinite(data[:, k]) * (error[:, k] != 0.0)
        D = np.sqrt(
            (xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid) ** 2
            + (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid) ** 2
        )
        W = np.exp(-0.5 / (seeing_fac / 2.35) ** 2 * D**2)
        zgrid[k, :, :] = (
            (data[sel, k][:, np.newaxis, np.newaxis] * W[sel]).sum(axis=0)
            / W[sel].sum(axis=0)
            * (scale**2 / area)
        )
    wi = np.searchsorted(wave, wstart, side="left")
    we = np.searchsorted(wave, wend, side="right")
    zimage = np.median(zgrid[wi : we + 1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


def get_objects(basefiles, attrs, full=False):
    """Extract selected FITS header attributes for a list of base files.

    Opens each FITS file within its tar archive, reads requested header
    attributes, and optionally appends guider-derived mirror illumination and
    throughput estimates.

    Args:
        basefiles (list[str]): List of FITS file paths inside tar archives.
        attrs (list[str]): FITS header keys to extract from the primary HDU.
        full (bool, optional): If True, append mirror illumination area and
            throughput per file. Defaults to False.

    Returns:
        list[list]: Nested list where each inner list contains the requested
            attributes (and optionally area, throughput) for the corresponding
            file in ``basefiles``.
    """
    s = []
    for fn in basefiles:
        tarname = get_tarname_from_filename(fn)
        t = tarfile.open(tarname, "r")
        F = fits.open(t.extractfile("/".join(fn.split("/")[-4:])))
        s.append([])
        for att in attrs:
            s[-1].append(F[0].header[att])
        if full:
            from .routine import get_mirror_illumination_guider, get_throughput

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
    """Select a representative 5-element subset from a list.

    Returns the first, three approximately equally spaced middle elements,
    and the last element. If the list length is <= 5, returns the list as-is.

    Args:
        lst (list): Input list of any type.

    Returns:
        list: Representative subset with up to 5 elements.
    """
    if len(lst) <= 5:
        return lst
    step = (len(lst) - 1) / 4
    indices = [0, round(step), round(2 * step), round(3 * step), len(lst) - 1]
    return [lst[i] for i in indices]


def get_previous_night(daten):
    """Compute the previous night's date string.

    Args:
        daten (str): Current date in YYYYMMDD format.

    Returns:
        str: Previous date in YYYYMMDD format.
    """
    daten_ = datetime(int(daten[:4]), int(daten[4:6]), int(daten[6:])) - timedelta(
        days=1
    )
    return f"{daten_.year:04d}{daten_.month:02d}{daten_.day:02d}"


def get_script_path():
    """Return the absolute directory of the running script.

    Uses sys.argv[0] resolved to an absolute path. This mirrors the legacy
    helper and is useful when scripts need to locate resources relative to
    their own file.

    Returns:
        str: Absolute path to the directory containing the current script.
    """
    return op.dirname(op.realpath(sys.argv[0]))


def count_matches(lines, loc, fib, cnt=5):
    """Score edge-anchored linear mappings to count peak/line matches.

    Emulates the legacy heuristic that tries simple linear mappings between the
    detected peak locations in a seed fiber and the reference line list by
    anchoring the first/last few detected peaks to the first/last reference
    lines. For each combination, counts how many lines can be matched within a
    small tolerance and returns the indices yielding the maximum count.

    Args:
        lines (astropy.table.Table): Reference line list with column 'col2'
            giving expected pixel/column positions.
        loc (list[np.ndarray]): Detected peak locations per fiber; only
            ``loc[fib]`` is used.
        fib (int): Fiber index to evaluate.
        cnt (int, optional): Number of candidate edges to try from each end.
            Defaults to 5.

    Returns:
        tuple[int, int]: Indices (k, j) maximizing the match count, where k is
        the offset from the start and j from the end in ``loc[fib]`` used to
        anchor the mapping.
    """
    x = lines["col2"]
    M = np.zeros((cnt, cnt))
    peaks = np.array(loc[fib])
    if peaks.size == 0 or len(x) == 0:
        return (0, 0)
    for k in np.arange(cnt):
        for j in np.arange(cnt):
            i1 = 0 + int(k)
            i2 = -1 - int(j)
            if i1 >= len(peaks) or abs(i2) > len(peaks):
                continue
            diff0 = [peaks[i1] - x[0], peaks[i2] - x[-1]]
            m = (diff0[1] - diff0[0]) / (x[-1] - x[0]) if (x[-1] - x[0]) != 0 else 0.0
            y = m * (x - x[0]) + diff0[0] + x
            count = 0
            for col in y:
                d = np.abs(col - peaks)
                if d.min() < 3.0:
                    count += 1
            M[k, j] = count
    return tuple(np.unravel_index(np.argmax(M), M.shape))


def get_standard_star_params(data, commonwave, xloc, yloc):
    """Estimate centroid, size, and DAR trends from a standard-star cube.

    Splits the wavelength range into 11 chunks, computes median fiber flux per
    chunk, fits a 2D Gaussian to the brightest region to estimate centroid and
    size per chunk, and then fits a quadratic trend of centroid vs. wavelength
    to derive differential atmospheric refraction (DAR) offsets across the full
    wavelength grid.

    Args:
        data: 2D array (Nfibers, Nwave) of rectified spectra for a standard star.
        commonwave: 1D wavelength grid corresponding to columns of ``data``.
        xloc: 1D array of fiber x positions.
        yloc: 1D array of fiber y positions.

    Returns:
        tuple: (x_center, y_center, xstd, ystd, xoff, yoff), where
        - x_center, y_center: float centroids at the central wavelength,
        - xstd, ystd: 1D arrays of Gaussian sigma along x and y (constant-valued),
        - xoff, yoff: 1D arrays of DAR offsets relative to the center column.
    """
    G = Gaussian2D()
    fitter = LevMarLSQFitter()
    wchunk = np.array([np.mean(chunk) for chunk in np.array_split(commonwave, 11)])
    dchunk = [np.median(chunk, axis=1) for chunk in np.array_split(data, 11, axis=1)]
    xc = 0.0 * wchunk
    yc = 0.0 * wchunk
    xs = 0.0 * wchunk
    ys = 0.0 * wchunk
    for i in np.arange(11):
        y = dchunk[i]
        ind = int(np.argmax(y))
        dist = np.sqrt((xloc - xloc[ind]) ** 2 + (yloc - yloc[ind]) ** 2)
        inds = dist < 3.0
        x_centroid = float(np.sum(y[inds] * xloc[inds]) / np.sum(y[inds]))
        y_centroid = float(np.sum(y[inds] * yloc[inds]) / np.sum(y[inds]))
        G.amplitude.value = float(y[ind])
        G.x_mean.value = x_centroid
        G.y_mean.value = y_centroid
        fit = fitter(G, xloc[inds], yloc[inds], y[inds])
        xc[i] = fit.x_mean.value
        yc[i] = fit.y_mean.value
        xs[i] = fit.x_stddev.value
        ys[i] = fit.y_stddev.value
    sel = xs > 0.0
    xoff = np.polyval(np.polyfit(wchunk[sel], xc[sel], 2), commonwave)
    yoff = np.polyval(np.polyfit(wchunk[sel], yc[sel], 2), commonwave)
    xstd = np.median(xs[sel]) * np.ones(commonwave.shape)
    ystd = np.median(ys[sel]) * np.ones(commonwave.shape)
    N = len(commonwave)
    mid = N // 2
    return xoff[mid], yoff[mid], xstd, ystd, xoff - xoff[mid], yoff - yoff[mid]


def get_bigarray(xloc, yloc):
    """Build a tiled grid of IFU fiber positions expanded in all directions.

    Starting from the observed fiber positions (xloc, yloc), extend the grid by
    repeating the edge rows and columns outward by roughly the native fiber
    spacing to create a larger lattice. Used for PSF normalization over a
    larger aperture when extracting a compact source.

    Args:
        xloc (np.ndarray): 1D x positions of fibers.
        yloc (np.ndarray): 1D y positions of fibers.

    Returns:
        tuple[np.ndarray, np.ndarray]: Arrays (BigX, BigY) of the expanded grid
        coordinates.
    """
    BigX = [xloc]
    BigY = [yloc]
    uy = np.unique(yloc)
    if len(uy) > 1:
        dy = float(np.mean(np.diff(uy)))
    else:
        dy = 0.0
    for _ in np.arange(1, 10):
        ny = dy + np.max(np.hstack(BigY))
        x = (
            np.hstack(BigX)[np.where(uy[-2] == np.hstack(BigY))[0]]
            if len(uy) > 1
            else np.hstack(BigX)
        )
        y = ny * np.ones(x.shape)
        BigY.append(y)
        BigX.append(x)
        uy = np.unique(np.hstack(BigY))
        ny = -dy + np.min(np.hstack(BigY))
        x = (
            np.hstack(BigX)[np.where(uy[1] == np.hstack(BigY))[0]]
            if len(uy) > 1
            else np.hstack(BigX)
        )
        y = ny * np.ones(x.shape)
        BigY.append(y)
        BigX.append(x)
        uy = np.unique(np.hstack(BigY))
    BigX = np.hstack(BigX)
    BigY = np.hstack(BigY)
    uy = np.unique(BigY)
    NX, NY = ([BigX], [BigY])
    for i in uy:
        sel = np.where(i == BigY)[0]
        if len(sel) < 2:
            continue
        dx = float(np.abs(np.mean(np.diff(BigX[sel]))))
        xn = np.min(BigX[sel]) - np.arange(1, 10) * dx
        xn2 = np.max(BigX[sel]) + np.arange(1, 10) * dx
        yn = i * np.ones(xn.shape)
        yn2 = i * np.ones(xn2.shape)
        NX.append(xn)
        NX.append(xn2)
        NY.append(yn)
        NY.append(yn2)
    BigX = np.hstack(NX)
    BigY = np.hstack(NY)
    return BigX, BigY


def robust_polyfit(x, y, order=3, niter=3):
    """MAD-clipped polynomial fit evaluated on x.

    Args:
        x: 1D array of x-coordinates.
        y: 1D array of y-values.
        order: Polynomial order.
        niter: Number of sigma-clipping iterations.

    Returns:
        Fitted y model evaluated at x.
    """
    sel = y > 0.0
    ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    for _ in np.arange(niter):
        a = np.abs(y - ymod)
        mad = np.median(a)
        sel = a < 3.0 * mad
        if sel.sum() > (order + 2):
            ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    return ymod
