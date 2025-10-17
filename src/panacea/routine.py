"""General routines orchestrating reduction steps.

Functions migrated per function_map.md.
"""
from __future__ import annotations

import os
import os.path as op
import glob
import tarfile
from datetime import datetime
import warnings
import logging

import numpy as np
from typing import List, Tuple, Optional
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import interp1d, interp2d
from scipy.signal import savgol_filter
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter

from .io import get_tarname_from_filename, get_filenames_from_tarfolder
from .ccd import base_reduction, get_powerlaw
from .fiber import get_spectra, weighted_extraction, modify_spectrum
from .trace import get_trace_shift
from .utils import check_if_standard, truncate_list
from .utils import build_weight_matrix, mask_skylines_cosmics


def get_ifucenfile(side: str, amp: str, virusconfig: str = "/work/03946/hetdex/maverick/virus_config", skiprows: int = 4) -> np.ndarray:
    """Load IFU center positions for a given side and amplifier.

    Args:
        side: One of {"uv", "orange", "red", "farred"}.
        amp: Amplifier identifier ("LL", "LU", "RL", "RU").
        virusconfig: Base path to the VIRUS configuration directory.
        skiprows: Number of header rows to skip in the mapping text file.

    Returns:
        Array of (x, y) fiber positions for the requested amplifier half.
    """
    file_dict = {
        "uv": "LRS2_B_UV_mapping.txt",
        "orange": "LRS2_B_OR_mapping.txt",
        "red": "LRS2_R_NR_mapping.txt",
        "farred": "LRS2_R_FR_mapping.txt",
    }

    ifucen = np.loadtxt(op.join(virusconfig, "IFUcen_files", file_dict[side]), usecols=[0, 1, 2], skiprows=skiprows)

    if amp == "LL":
        return ifucen[140:, 1:3][::-1, :]
    if amp == "LU":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RL":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RU":
        return ifucen[140:, 1:3][::-1, :]
    raise ValueError(f"Unknown amp: {amp}")



def get_mirror_illumination(fn=None, default=51.4e4):
    """Compute mirror illumination area using hetillum (external), fallback to default.

    Args:
        fn: FITS filename or file-like to read header geometry parameters.
        default: Default area (cm^2) to return on failure.

    Returns:
        Mirror illumination area in cm^2.
    """
    try:
        F = fits.open(fn)
        names = ['RHO_STRT', 'THE_STRT', 'PHI_STRT', 'X_STRT', 'Y_STRT']
        r, t, p, x, y = [F[0].header[name] for name in names]
        cmd = (
            '/work/03730/gregz/maverick/illum_lib/hetillum -p'
            f' -x "[{x:.4f},{y:.4f},{p:.4f}]" "[{0.042:.4f},{0.014:.4f}]" 256'
        )
        mirror_illum = float(os.popen(cmd).read().split('\n')[0])
        area = mirror_illum * default
    except Exception:
        area = default
    return area


def get_mirror_illumination_guider(fn, exptime, default=51.4e4, path='/work/03946/hetdex/maverick'):
    """Estimate mirror illumination using guider tarballs around exposure time.

    Args:
        fn: Science filename to infer datetime.
        exptime: Exposure time (s).
        default: Default area if guider data not found.
        path: Base path containing guider tarballs by date.
    Returns:
        area: Mirror illumination area in cm^2.
    """
    try:
        M = []
        # infer date directory
        f = op.basename(fn)
        DT = f.split('_')[0]
        y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
        d0 = datetime(y, m, d, h, mi, s)
        tarfolder = op.join(path, 'gc1', '*.tar')
        tarfolder = glob.glob(tarfolder)
        if len(tarfolder) == 0:
            return default
        T = tarfile.open(tarfolder[0], 'r')
        init_list = sorted([name for name in T.getnames() if name.endswith('.fits')])
        final_list = []
        for t in init_list:
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.0) * (p < exptime + 10.0):
                final_list.append(t)
        # truncate list to 5
        final_list = truncate_list(final_list)
        for f2 in final_list:
            fobj = T.extractfile(T.getmember(f2))
            M.append(get_mirror_illumination(fobj))
        M = np.array(M)
        sel = M != default
        if sel.sum() > 0.0:
            area = float(np.mean(M[sel]))
        else:
            area = default
        return area
    except Exception:
        return default


def get_throughput(fn, exptime, path='/work/03946/hetdex/maverick'):
    """Estimate atmospheric throughput from guider telemetry.

    Computes the average of the TRANSPAR header keyword from guider cameras
    GC1 and GC2 over the time window that overlaps the science exposure. Only
    frames with GUIDLOOP == 'ACTIVE' are considered. If no valid frames are
    found, or if the resulting average is out of a reasonable range, the
    function returns 1.0 as a safe default.

    Args:
        fn (str): Science filename used to infer the mid/exposure start time.
            The timestamp is parsed from the basename as YYYYMMDD_hhmmss.
        exptime (float): Exposure time in seconds; defines the time window
            for selecting guider frames (from t0-10s to t0+exptime+10s).
        path (str, optional): Base directory containing guider tarballs,
            expected at "{path}/gc1/gc1.tar" and "{path}/gc2/gc2.tar".
            Defaults to '/work/03946/hetdex/maverick'.

    Returns:
        float: Mean TRANSPAR value in [0.1, 1.1] from active guider frames. If
        the computed value is NaN, outside the allowed range, or no frames are
        available, returns 1.0.
    """
    attr = ['GUIDLOOP', 'MJD', 'TRANSPAR']
    M = []
    f = op.basename(fn)
    DT = f.split('_')[0]
    y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
    d0 = datetime(y, m, d, h, mi, s)
    for gp in ['gc1', 'gc2']:
        tarfolder = op.join(path, gp, f'{gp}.tar')
        T = tarfile.open(tarfolder, 'r')
        init_list = sorted([name for name in T.getnames() if name.endswith('.fits')])
        final_list = []
        for t in init_list:
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.0) * (p < exptime + 10.0):
                final_list.append(t)
        final_list = truncate_list(final_list)
        for fnt in final_list:
            fobj = T.extractfile(T.getmember(fnt))
            f1 = fits.open(fobj)
            if f1[1].header.get('GUIDLOOP') == 'ACTIVE':
                M.append([])
                for att in attr:
                    M[-1].append(f1[1].header.get(att, 0))
    throughput = np.zeros((len(M),)) if len(M) else np.zeros((1,))
    for i, mi in enumerate(M):
        if len(mi) < 2:
            continue
        if mi[2] > 0.0:
            throughput[i] = mi[2]
    t = float(np.mean(throughput[throughput > 0.0])) if throughput.size else 1.0
    if not np.isfinite(t) or t > 1.1 or t < 0.1:
        t = 1.0
    return t


def extract_sci(sci_path: str, amps: List[str], flat: np.ndarray, array_trace: np.ndarray,
                array_wave: np.ndarray, bigW: np.ndarray, masterbias: np.ndarray,
                pos: np.ndarray, commonwave: np.ndarray):
    """Extract, sky-normalize, and rectify science spectra for a two-amp exposure.

    This function mirrors the legacy extract_sci workflow while explicitly
    requiring a common wavelength grid for rectification. For each exposure in
    the pair of amplifier halves, it performs basic CCD reduction, subtracts the
    master bias, estimates and removes a low-order power-law background, and
    carries out a weighted extraction using the provided flat and trace models.
    Spectra are then adjusted via ``modify_spectrum`` and linearly interpolated
    onto ``commonwave`` with propagated errors.

    Args:
        sci_path: Path pattern pointing to the LL amplifier file inside the tar; the
            corresponding other amp in ``amps`` will be substituted for each exposure.
        amps: Two amplifier labels corresponding to the pair to combine (e.g.,
            ["LL", "LU"] or ["RL", "RU"]). The LL in ``sci_path`` is replaced by
            each of these to collect files.
        flat: 2D flat-field image used during weighted extraction.
        array_trace: 2D array of trace center positions per fiber and column.
        array_wave: 2D wavelength solution array per fiber and column.
        bigW: Unused placeholder retained for API compatibility with legacy code.
        masterbias: 2D master bias frame to subtract from raw images.
        pos: 2D array of fiber positions; columns [x, y] are used.
        commonwave: 1D wavelength grid onto which spectra are rectified.

    Returns:
        Tuple containing:
        - images (np.ndarray): Stack of bias-subtracted 2D images used in extraction.
        - spec_list (np.ndarray): Rectified spectra on ``commonwave`` per exposure
          with shape (Nexp, Nfibers, Nwave).
        - orig_list (np.ndarray): Pre-rectification spectra on native wavelength.
        - clist (np.ndarray): Continuum normalization coefficients from extraction.
        - flist (np.ndarray): Per-fiber flux list or normalization factors.
        - Flist (np.ndarray): Per-exposure 2D extraction weight images.
        - error_list (np.ndarray): Rectified 1-sigma errors aligned with ``spec_list``.
        - hdr_list (list[astropy.io.fits.Header]): FITS headers from the reduced files.
    """
    log = logging.getLogger(__name__)
    files1 = get_filenames_from_tarfolder(get_tarname_from_filename(sci_path), sci_path.replace('LL', amps[0]))
    files2 = get_filenames_from_tarfolder(get_tarname_from_filename(sci_path), sci_path.replace('LL', amps[1]))
    tarnames = [get_tarname_from_filename(file) for file in files1]
    xloc, yloc = (pos[:, 0], pos[:, 1])
    array_list, hdr_list = ([], [])
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        log.info('Prepping sci %s', filename1)
        array_flt1, e1, header = base_reduction(filename1, tarname=tarname, get_header=True)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
        hdr_list.append(header)
    sci_array = np.sum(array_list, axis=0) if len(array_list) > 1 else np.squeeze(np.array(array_list))
    Xx = np.arange(flat.shape[1])
    Yx = np.arange(flat.shape[0])
    I2d = interp2d(Xx, Yx, flat, kind='cubic', bounds_error=False, fill_value=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        shifts = get_trace_shift(sci_array, flat, array_trace, Yx)
    flat = I2d(Xx, Yx + shifts)
    log.info('Found trace shift median: %0.3f', float(np.median(shifts)))
    spec_list, error_list, orig_list = ([], [], [])
    clist, flist, Flist = ([], [], [])
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        log.info('Fiber extraction sci %s', filename1)
        array_flt1, e1 = base_reduction(filename1, tarname=tarname)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_err = np.vstack([e1, e2])
        array_flt[:] -= masterbias
        spectrum = get_spectra(array_flt, array_trace)
        plaw, norm = get_powerlaw(array_flt, array_trace, spectrum)
        array_flt[:] -= plaw
        spectrum, error, c, fl, Fimage = weighted_extraction(array_flt, array_err, flat, array_trace)
        sel = ~np.isfinite(spectrum)
        spectrum[sel] = 0.0
        error[sel] = 0.0
        spectrum, error = modify_spectrum(spectrum, error, array_wave)
        # rectify to commonwave
        speclist, errorlist = ([], [])
        for fiber in np.arange(array_wave.shape[0]):
            I = interp1d(array_wave[fiber], spectrum[fiber], kind='linear', fill_value='extrapolate')
            # propagate error via coefficients of linear interp on commonwave grid
            coV = np.zeros((len(commonwave), spectrum.shape[1]))
            for i in np.arange(len(commonwave)):
                w = np.zeros((len(commonwave),))
                w[i] = 1.0
                coV[i] = np.interp(array_wave[fiber], commonwave, w)
            error_interp = np.sqrt((coV * error[fiber] ** 2).sum(axis=1))
            sel = error[fiber] == 0.0
            if sel.sum() > 0.0:
                nsel = coV[:, sel].sum(axis=1) > 0.0
                error_interp[nsel] = 0.0
            speclist.append(I(commonwave))
            errorlist.append(error_interp)
        spec_list.append(np.array(speclist))
        error_list.append(np.array(errorlist))
        orig_list.append(spectrum)
        clist.append(c)
        flist.append(fl)
        Flist.append(Fimage)
    images = np.array(array_list)
    return images, np.array(spec_list), np.array(orig_list), np.array(clist, dtype=float), np.array(flist), np.array(Flist), np.array(error_list), hdr_list


def find_source(dx: np.ndarray, dy: np.ndarray, skysub: np.ndarray, commonwave: np.ndarray,
                obj: str, specn: str, error: np.ndarray, xoff: np.ndarray, yoff: np.ndarray,
                wave_0: float, ispec: np.ndarray):
    """Locate a compact source using spatial and spectral convolution with S/N selection.

    This routine applies a Gaussian convolution along wavelength and a
    Gaussian-weighted sum across neighboring fibers after masking skylines and
    likely cosmic rays to form a convolved signal and error image. It then
    searches for the wavelength bin yielding the highest median S/N across
    fibers, collapses a narrow window around this location, and estimates the
    source centroid and seeing by fitting a 2D Gaussian to nearby fibers. A
    simple S/N threshold is used to accept or reject a detection.

    Args:
        dx: 1D array of fiber x positions.
        dy: 1D array of fiber y positions.
        skysub: 2D sky-subtracted rectified spectra (Nfibers, Nwave).
        commonwave: 1D common wavelength grid corresponding to columns of ``skysub``.
        obj: Object name (used only for logging).
        specn: Spectral channel name (e.g., 'uv', 'orange', 'red', 'farred'); used for masking.
        error: 2D error array aligned with ``skysub``.
        xoff: 1D DAR x offsets per wavelength index; will be recentered to the detection column.
        yoff: 1D DAR y offsets per wavelength index; will be recentered to the detection column.
        wave_0: Central wavelength reference (unused here; retained for API compatibility).
        ispec: 2D intermediate spectrum used to detect cosmics (same shape as ``skysub``).

    Returns:
        None if no significant source is found (S/N <= 5). Otherwise a tuple:
        - x_centroid (float): Source x centroid in IFU coordinates.
        - y_centroid (float): Source y centroid in IFU coordinates.
        - x_std (np.ndarray): 1D array of Gaussian sigma along x per wavelength (constant-valued).
        - y_std (np.ndarray): 1D array of Gaussian sigma along y per wavelength (constant-valued).
        - xoff (np.ndarray): Recentered DAR x offsets with zero at detection column.
        - yoff (np.ndarray): Recentered DAR y offsets with zero at detection column.
    """


    log = logging.getLogger(__name__)
    D = np.sqrt((dx - dx[:, np.newaxis]) ** 2 + (dy - dy[:, np.newaxis]) ** 2)
    # Convolve spatially and along wavelength similar to legacy
    # Build weights for near-neighbor fibers
    W = build_weight_matrix(dx, dy, sig=0.75)
    # Mask skylines and cosmics
    mask = mask_skylines_cosmics(commonwave, skysub, specn, error)
    Z = skysub * 1.0
    E = error ** 2
    Z[mask] = np.nan
    E[mask] = np.nan
    ispec = ispec * 1.0
    T = ispec * 1.0
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], (D < 1.5).astype(float))
    YY = ispec / T
    YY[np.isnan(YY)] = 0.0
    Z[YY > 0.2] = np.nan
    E[YY > 0.2] = np.nan
    G = Gaussian1DKernel(1.5)
    for i in np.arange(skysub.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment='fill', fill_value=0.0)
        E[i, :] = convolve(E[i, :], G, nan_treatment='fill', fill_value=0.0)
    sderror = np.sqrt(E)
    Y = Z * 0.0
    sel = E > 0.0
    Y[sel] = Z[sel] / np.sqrt(E[sel])
    Y[~np.isfinite(Y)] = 0.0
    ind = np.unravel_index(np.nanargmax(Y[:, 50:-50], axis=None), Y[:, 50:-50].shape)
    loc = ind[1] + 50
    BN = 25
    dimage = np.sum(Z[:, (loc-BN):(loc+BN+1)], axis=1)
    derror = np.sqrt(np.sum(sderror[:, (loc-BN):(loc+BN+1)] ** 2, axis=1))
    sn = dimage * 0.0
    for i in np.arange(len(dimage)):
        sel = D[i, :] < 1.5
        S = np.sum(dimage[sel])
        N = np.sqrt(np.sum(derror[sel] ** 2))
        sn[i] = S / N if N > 0 else 0.0
    SN = float(np.nanmax(sn)) if np.isfinite(sn).any() else 0.0
    if SN <= 5.0:
        log.info('%s, %s: No source found, s/n too low: %0.2f', obj, specn, SN)
        return None
    ind = int(np.argmax(dimage))
    dist = np.sqrt((dx - dx[ind]) ** 2 + (dy - dy[ind]) ** 2)
    inds = dist < 1.5
    x_centroid = float(np.sum(dimage[inds] * dx[inds]) / np.sum(dimage[inds]))
    y_centroid = float(np.sum(dimage[inds] * dy[inds]) / np.sum(dimage[inds]))
    G2 = Gaussian2D()
    fitter = LevMarLSQFitter()
    G2.amplitude.value = float(dimage[ind])
    G2.x_mean.value = x_centroid
    G2.y_mean.value = y_centroid
    fit = fitter(G2, dx[inds], dy[inds], dimage[inds])
    seeing = 2.35 * np.sqrt(fit.x_stddev * fit.y_stddev)
    if seeing < 0.75:
        log.info('%s, %s: source s/n %0.2f rejected for too small FWHM, col: %i', obj, specn, SN, loc)
        return None
    log.info('%s, %s: source found at s/n: %0.2f, fwhm: %s', obj, specn, SN, seeing)
    X = np.ones(commonwave.shape)
    xoff = xoff - xoff[loc]
    yoff = yoff - yoff[loc]
    return x_centroid, y_centroid, fit.x_stddev.value * X, fit.y_stddev.value * X, xoff, yoff



def fit_response_cont(wv: np.ndarray, sky: np.ndarray, skip: int = 5, fil_len: int = 95, func=np.array):
    """Estimate a smooth continuum for response calibration.

    Applies a Savitzky–Golay filter to the input sky spectrum with iterative
    outlier masking around sharp spectral features. Intended to approximate the
    broadband continuum for computing response functions.

    Args:
        wv (np.ndarray): 1D wavelength vector.
        sky (np.ndarray): 1D sky spectrum sampled at ``wv``.
        skip (int, optional): Number of pixels to expand when masking around
            detected sharp features. Defaults to 5.
        fil_len (int, optional): Window length for Savitzky–Golay filter. Must
            be odd and > polyorder. Defaults to 95.
        func (callable, optional): Function applied to residuals to detect
            negative outliers; pass ``np.array`` for identity. Defaults to np.array.

    Returns:
        np.ndarray: Smoothed continuum estimate sampled at ``wv``.
    """
    skym_s = 1.0 * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    y = np.abs(np.diff(np.hstack([sky[0], sky])))
    sel = np.where(y > 1.5 * np.median(sky))[0]
    for j in np.arange(1, skip + 1):
        sel = np.union1d(sel, sel + 1)
        sel = np.union1d(sel, sel - 1)
    sel = np.sort(np.unique(sel))
    sel = sel[skip:-skip]
    good = np.setdiff1d(allind, sel)
    skym_s = 1.0 * sky
    skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    for _ in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        outlier = func(sky - sky_sm) < -1.5 * mad
        sel = np.where(outlier)[0]
        for j in np.arange(1, skip + 1):
            sel = np.union1d(sel, sel + 1)
            sel = np.union1d(sel, sel - 1)
        sel = np.sort(np.unique(sel))
        sel = sel[skip:-skip]
        good = np.setdiff1d(allind, sel)
        skym_s = 1.0 * sky
        skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    return sky_sm


def get_response(objname: str, commonwave: np.ndarray, spec: np.ndarray, specname: str):
    """Compute a scalar response vector from a spectrophotometric standard.

    If the target name matches a known standard star, this function loads the
    corresponding reference AB magnitudes, converts them to f_lambda, and
    derives a smooth continuum ratio between the observed spectrum and the
    model via ``fit_response_cont``. The returned array is a multiplicative
    response to bring the observed spectrum onto an absolute scale. If the
    object is not recognized as a standard or the file is missing, returns None.

    Args:
        objname: Object name to match against known standards (case-insensitive,
            substring match).
        commonwave: 1D wavelength grid of the observed ``spec``.
        spec: 1D observed spectrum of the standard sampled at ``commonwave``.
        specname: Spectral channel name (unused; retained for logging/compatibility).

    Returns:
        np.ndarray | None: Multiplicative response array sampled at ``commonwave``
        if a standard is identified; otherwise None.
    """
    if not check_if_standard(objname):
        return None
    # External standards directory; keep path for compatibility
    filename = op.join('/work/03946/hetdex/maverick/virus_config/standards', f'm{objname.lower()}.dat.txt')
    try:
        wave, standardmag = np.loadtxt(filename, usecols=(0, 1), unpack=True)
    except Exception:
        return None
    fnu = 10 ** (0.4 * (-48.6 - standardmag))
    standard_flam = fnu * 2.99792e18 / wave ** 2
    standard_wave = wave
    flam = np.interp(commonwave, standard_wave, standard_flam)
    cont = fit_response_cont(commonwave, spec / flam, fil_len=11)
    return 1.0 / cont
