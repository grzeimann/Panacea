"""CCD-level processing functions for Panacea.

Functions copied from full_lrs2_reduction.py as per function_map.md.
Docstrings are updated to Google style where applicable.
"""

import warnings
import tarfile
import logging

import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location
from scipy.interpolate import interp1d, griddata
from scipy.signal import savgol_filter

from .utils import power_law
from .io import get_tarname_from_filename
from .sky import make_avg_spec

log = logging.getLogger(__name__)


def orient_image(image, amp, ampname):
    """Adjust the orientation of an image array based on amplifier info.

    The transformation is performed in-place and also returned for convenience.

    Args:
        image: 2D image array to be re-oriented.
        amp: Amplifier identifier (e.g., "LU", "RL").
        ampname: Optional AMPNAME keyword (e.g., "LR", "UL").

    Returns:
        The re-oriented image array (same object as input, modified in-place).
    """
    if amp == "LU":
        image[:] = image[::-1, ::-1]
    if amp == "RL":
        image[:] = image[::-1, ::-1]
    if ampname is not None:
        if ampname == "LR" or ampname == "UL":
            image[:] = image[:, ::-1]
    return image


def find_cosmics(Y, E, trace, thresh=8.0, **kwargs):
    """Detect cosmic ray hits near fiber traces using deviation thresholding.

    Args:
        Y: 2D image values.
        E: 2D error/variance estimates aligned with ``Y``.
        trace: Trace positions (fibers x columns).
        thresh: Sigma-like threshold; larger values are less sensitive.

    Returns:
        Boolean mask with True indicating suspected cosmic ray pixels.
    """
    x = np.arange(trace.shape[1])
    C = Y * 0.0
    for fiber in np.arange(trace.shape[0]):
        indl = np.floor(trace[fiber]).astype(int)
        T = np.zeros((4, trace.shape[1], 4))
        flag = True
        for ss, k in enumerate(np.arange(-1, 3)):
            try:
                T[0, :, ss] = Y[indl + k, x]
                T[1, :, ss] = E[indl + k, x]
                T[2, :, ss] = indl + k
                T[3, :, ss] = x
            except Exception:
                flag = False
        if flag:
            m = np.median(T[0], axis=1)
            P = np.abs(T[0] - m[:, np.newaxis]) / T[1]
            C[T[2][P > thresh].astype(int), T[3][P > thresh].astype(int)] = 1.0
    C = np.array(C, dtype=bool)
    log.info("Number of fiber pixels hit by cosmics: %i", C.sum())
    return C


def base_reduction(filename, tarname=None, get_header=False):
    """Basic CCD reduction: overscan subtraction, trim, gain, noise.

    Reads a FITS file from disk or from within a TAR archive, applies overscan
    subtraction, trims the overscan columns, applies gain, and computes a simple
    noise estimate.

    Args:
        filename: Path to FITS file (or path inside archive when ``tarname`` given).
        tarname: Optional path to a TAR archive containing the file.
        get_header: If True, also return the FITS header object.

    Returns:
        If ``get_header`` is False: (image, error)
        If ``get_header`` is True: (image, error, header)

    Notes:
        In case of failure to open an archived file, returns zero arrays with the
        expected CCD shape (1032, 2064).
    """
    if tarname is None:
        a = fits.open(filename)
    else:
        try:
            t = tarfile.open(tarname, "r")
            a = fits.open(t.extractfile("/".join(filename.split("/")[-4:])))
        except Exception:
            # Fallback to zeros if we cannot open the archived path
            return np.zeros((1032, 2064)), np.zeros((1032, 2064))
    image = np.array(a[0].data, dtype=float)
    # overscan subtraction
    overscan_length = int(32 * (image.shape[1] / 1064))
    O = biweight_location(image[:, -int(overscan_length - 2) :])
    image[:] = image - O
    # trim image
    image = image[:, :-overscan_length]
    gain = a[0].header.get("GAIN", 0.85)
    gain = np.where(gain > 0.0, gain, 0.85)
    rdnoise = a[0].header.get("RDNOISE", 3.0)
    rdnoise = np.where(rdnoise > 0.0, rdnoise, 3.0)
    amp = (a[0].header["CCDPOS"].replace(" ", "") + a[0].header["CCDHALF"].replace(" ", ""))
    ampname = a[0].header.get("AMPNAME", None)
    header = a[0].header
    aimg = orient_image(image, amp, ampname) * gain
    E = np.sqrt(rdnoise**2 + np.where(aimg > 0.0, aimg, 0.0))
    if tarname is not None:
        t.close()
    if get_header:
        return aimg, E, header
    return aimg, E


def get_powerlaw(image, trace, spec):
    """Calculate power-law normalization model near fiber traces.

    Args:
        image: 2D CCD image.
        trace: Fiber trace positions (fibers x columns).
        spec: Per-fiber spectra aligned with ``trace``.

    Returns:
        Tuple of (model_image, column_normalization).
    """
    YM, XM = np.indices(image.shape)
    inds = []
    for j in np.arange(trace.shape[0]):
        inds.append(np.where(np.abs(YM - trace[j, np.newaxis, :]).ravel() < 7.0)[0])
    inds = np.hstack(inds)
    inds = np.unique(inds)
    inds = np.setdiff1d(np.arange(image.shape[0] * image.shape[1]), inds)
    y, x = np.unravel_index(inds, image.shape)
    xlim = [0, image.shape[1]]
    xy = np.nanmedian(spec, axis=1)
    bottom = np.abs(np.nanpercentile(xy, 15))
    sel = np.where(xy > (15.0 * bottom))[0]
    xp = np.hstack([np.arange(xlim[0], xlim[1], 128), image.shape[1] - 1])
    ylim = [0, image.shape[0]]
    yz = np.hstack([np.arange(ylim[0], ylim[1], 128), image.shape[0] - 1])
    plaw, XX, YY = ([], [], [])
    YM2, XM2 = np.indices(trace.shape)
    for xi in xp:
        for yi in yz:
            d = np.sqrt((yi - trace) ** 2 + (xi - XM2) ** 2)
            plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2.0, c4=1.0, sig=1.5)))
            XX.append(xi)
            YY.append(yi)
    for xi in xp:
        for s in sel:
            y0 = int(trace[s, xi])
            for i in np.arange(-4, 6, 2):
                d = np.sqrt(((y0 + i) - trace) ** 2 + (xi - XM2) ** 2)
                plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2.0, c4=1.0, sig=1.5)))
                YY.append(y0 + i)
                XX.append(xi)
    plaw, XX, YY = [np.hstack(j) for j in [plaw, XX, YY]]
    grid_x, grid_y = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
    C = griddata(np.array([XX, YY]).swapaxes(0, 1), plaw, (grid_x, grid_y), method="cubic", fill_value=0.0)
    norm = np.zeros((image.shape[1],))
    for b in np.arange(image.shape[1]):
        selb = x == b
        norm[b] = np.median(image[y[selb], x[selb]] / C[y[selb], x[selb]])
    return C * savgol_filter(norm, 141, 1)[np.newaxis, :], norm


def get_bigW(amp, array_wave, array_trace, image):
    """Compute the wavelength grid per pixel over the full CCD.

    Args:
        amp: Unused (kept for compatibility with legacy signature).
        array_wave: Per-fiber wavelength solutions (fibers x columns).
        array_trace: Per-fiber trace positions.
        image: Reference 2D image providing the output shape.

    Returns:
        bigW image with wavelength mapped at each CCD pixel row for each column.
    """
    bigW = np.zeros(image.shape)
    Y, X = np.indices(array_wave.shape)
    YY, XX = np.indices(image.shape)
    for x, at, aw, xx, yy in zip(
        np.array_split(X, 2, axis=0),
        np.array_split(array_trace, 2, axis=0),
        np.array_split(array_wave, 2, axis=0),
        np.array_split(XX, 2, axis=0),
        np.array_split(YY, 2, axis=0),
    ):
        for j in np.arange(at.shape[1]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p0 = np.polyfit(at[:, j], aw[:, j], 7)
            bigW[yy[:, j], j] = np.polyval(p0, yy[:, j])
    return bigW


def get_bigF(array_trace, image):
    """Compute a smoothed spatial profile model across the CCD.

    Args:
        array_trace: Per-fiber trace positions.
        image: Reference image providing the output shape.

    Returns:
        2D model of the expected trace center locations per column.
    """
    bigF = np.zeros(image.shape)
    Y, X = np.indices(array_trace.shape)
    YY, XX = np.indices(image.shape)
    n, m = array_trace.shape
    F0 = array_trace[:, int(m / 2)]
    for j in np.arange(image.shape[1]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = np.polyfit(array_trace[:, j], F0, 7)
        bigF[:, j] = np.polyval(p0, YY[:, j])
    return bigF


def get_masterbias(files, amp):
    """Compute a master bias from a list of raw frames.

    Args:
        files: List of file templates containing 'LL' to be replaced by amp.
        amp: Amplifier identifier used to select files.

    Returns:
        2D master bias image.
    """
    files = [file.replace("LL", amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]

    biassum = np.zeros((len(files), 1032, 2064))
    for j, filename in enumerate(files):
        tarname = tarnames[j]
        a, error = base_reduction(filename, tarname=tarname)
        biassum[j] = a
    return biweight_location(biassum, axis=0)


def get_masterarc(files, amp, arc_names, masterbias, specname, trace):
    """Build a master arc image by combining suitable arc exposures.

    Args:
        files: List of file templates (will replace 'LL' with ``amp``).
        amp: Amplifier identifier.
        arc_names: Lowercase object names that indicate arc exposures.
        masterbias: Master bias image.
        specname: Channel name (kept for signature compatibility).
        trace: Trace information (currently unused here).

    Returns:
        2D master arc image.
    """
    files = [file.replace("LL", amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    arcsum = np.zeros((1032, 2064))
    cnt = np.zeros((1032, 2064))
    for filename, tarname in zip(files, tarnames):
        t = tarfile.open(tarname, "r")
        f = fits.open(t.extractfile("/".join(filename.split("/")[-4:])))
        if f[0].header["OBJECT"].lower() in arc_names:
            a, e = base_reduction(filename, tarname=tarname)
            a[:] -= masterbias
            if np.median(a) < 3000.0:
                arcsum += a
                cnt += 1.0
    return arcsum / cnt


def get_mastertwi(files, amp, masterbias):
    """Compute a master twilight flat image.

    Args:
        files: List of file templates (will replace 'LL' with ``amp``).
        amp: Amplifier identifier.
        masterbias: Master bias image.

    Returns:
        2D master twilight flat (median-normalized per exposure before combine).
    """
    files = [file.replace("LL", amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    listtwi = []
    for filename, tarname in zip(files, tarnames):
        a, e = base_reduction(filename, tarname=tarname)
        a[:] -= masterbias
        if np.percentile(a, 75) > 100.0:
            listtwi.append(a)
    twi_array = np.array(listtwi, dtype=float)
    norm = np.median(twi_array, axis=(1, 2))[:, np.newaxis, np.newaxis]
    return np.median(twi_array / norm, axis=0)


def get_twiflat_field(files, amps, array_wave, array_trace, bigW, masterbias, specname):
    """Create a normalized flat-field image from twilight/internal flats.

    This computes fiber spectra from flats, removes a smoothly varying power-law
    component, derives fiber-to-fiber (ftf) normalization, and returns a flat
    image normalized along fibers and across wavelength.

    Args:
        files: File templates for two amplifiers.
        amps: Two amplifier labels (e.g., ["LL", "LU"]).
        array_wave: Wavelength grid (fibers x columns).
        array_trace: Trace positions (fibers x columns).
        bigW: High-resolution per-pixel wavelength grid for the 2D image.
        masterbias: Master bias frame to subtract.
        specname: Channel name for logging.

    Returns:
        2D normalized flat image.
    """
    files1 = [file.replace("LL", amps[0]) for file in files]
    files2 = [file.replace("LL", amps[1]) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    array_list = []
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        array_flt1, e1 = base_reduction(filename1, tarname=tarname)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
    array_list = np.array(array_list)
    if len(array_list) > 1:
        norm = np.median(array_list, axis=(1, 2))
        array_flt = np.median(array_list / norm[:, np.newaxis, np.newaxis], axis=0)
    else:
        array_flt = np.squeeze(np.array(array_list))
        array_flt[:] /= np.median(array_flt)

    x = np.arange(array_wave.shape[1])
    spectrum = array_trace * 0.0
    for fiber in np.arange(array_wave.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        try:
            spectrum[fiber] = array_flt[indl, x] / 2.0 + array_flt[indh, x] / 2.0
        except Exception:
            spectrum[fiber] = 0.0

    plaw, norm = get_powerlaw(array_flt, array_trace, spectrum)
    array_flt[:] -= plaw
    array_flt[:] = np.where(array_flt < 0.0, 0.0, array_flt)

    smooth = savgol_filter(spectrum, 315, 1, axis=1)
    avg = biweight_location(smooth, axis=(0,))
    norm = biweight_location(smooth / avg, axis=(1,))
    nw, ns = make_avg_spec(array_wave, spectrum / norm[:, np.newaxis], binsize=41, per=50)
    I = interp1d(nw, ns, kind="linear", fill_value="extrapolate")
    ftf = spectrum * 0.0
    for fiber in np.arange(array_wave.shape[0]):
        model = I(array_wave[fiber])
        ftf[fiber] = savgol_filter(spectrum[fiber] / model, 151, 1)
    nw1, ns1 = make_avg_spec(array_wave, spectrum / ftf, binsize=41, per=50)

    I2 = interp1d(nw1, ns1, kind="quadratic", fill_value="extrapolate")
    modelimage = I2(bigW)
    flat = array_flt / modelimage
    flat[~np.isfinite(flat)] = 0.0
    flat[flat < 0.0] = 0.0
    return flat
