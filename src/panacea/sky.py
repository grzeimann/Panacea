"""Sky and statistical utilities.

Migrated functions per function_map.md.
"""

import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Polynomial2D
from .wavelength import find_peaks


def make_avg_spec(wave, spec, binsize=35, per=50):
    """Generate averaged spectrum by binning and taking percentiles.

    Args:
        wave: Wavelength array (any shape). Will be raveled for processing.
        spec: Spectrum array aligned with ``wave`` (same shape).
        binsize: Approximate number of samples per bin.
        per: Percentile to compute per bin (e.g., 50 for median).

    Returns:
        Tuple of (nwave, nspec) where both are 1D arrays of unique averaged
        wavelengths and their corresponding percentile flux values.
    """
    ind = np.argsort(wave.ravel())
    T = 1
    for p in wave.shape:
        T *= p
    wchunks = np.array_split(wave.ravel()[ind], T / binsize)
    schunks = np.array_split(spec.ravel()[ind], T / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    nspec = np.array([np.percentile(chunk, per) for chunk in schunks])
    nwave, nind = np.unique(nwave, return_index=True)
    return nwave, nspec[nind]


def sky_subtraction(rect, error, xloc, yloc):
    """Build a 2D sky model and per-column residual, returning sky+res.

    Algorithm mirrors legacy implementation: selects sky-dominated fibers,
    builds an initial median sky, detects sky-line regions, fits fiber scale
    factors via least-squares, then fits a 2D polynomial across IFU geometry
    to model spatial variation. Returns the modeled sky plus residual ramp.

    Args:
        rect: Rectified spectra (fibers x columns).
        error: Error array aligned with ``rect``.
        xloc: IFU x positions per fiber.
        yloc: IFU y positions per fiber.

    Returns:
        2D array of the sky model including residual ramp per column.
    """
    y = np.median(rect, axis=1)
    selg = y != 0.0
    v = np.percentile(y[selg], 5)
    init_sel = selg * (y < v)
    init = np.percentile(rect[init_sel], 50, axis=0)
    df = np.diff(init)
    df = np.hstack([df[0], df])
    cont = np.abs(df) < np.percentile(np.abs(df), 25)
    G = Gaussian1DKernel(15)
    tempy = init * 1.0
    tempy[~cont] = np.nan
    smooth_back = convolve(tempy, G, nan_treatment="interpolate", preserve_nan=False)
    peak_loc, sn, v = find_peaks(init - smooth_back, thresh=3)
    locs = np.round(peak_loc).astype(int)
    locs = np.sort(np.hstack([locs - 2, locs - 1, locs, locs + 1, locs + 2]))
    locs = locs[np.where(locs >= 0)[0]]
    locs = locs[np.where(locs < 2064)[0]]
    facs = np.arange(0.7, 1.3, 0.01)
    cnt = facs * 0.0
    sol = np.ones((280,))
    for k in np.arange(280):
        good = np.zeros(rect[k].shape, dtype=bool)
        good[locs] = True
        good[error[k] == 0.0] = False
        if good.sum() > 50.0:
            for j, i in enumerate(facs):
                cnt[j] = (((rect[k] - i * init) * init)[good]).sum()
            ind = np.argsort(cnt)
            sol[k] = np.interp(0.0, cnt[ind], facs[ind])
    sol[np.abs(sol - 1.3) < 0.01] = 1.0
    n1 = np.median(sol[:140])
    n2 = np.median(sol[140:])
    nsol = sol * 1.0
    nsol[:140] = sol[:140] / n1
    nsol[140:] = sol[140:] / n2
    fitter = LevMarLSQFitter()
    P = Polynomial2D(2)
    good = nsol != 0.0
    fit = fitter(P, xloc[good], yloc[good], nsol[good])
    off = np.abs(sol - fit(xloc, yloc))
    mad = np.median(off)
    good = (nsol != 0.0) * (off <= 2.0 * mad)
    fit = fitter(P, xloc[good], yloc[good], nsol[good])
    model = fit(xloc, yloc)
    model[:140] *= n1
    model[140:] *= n2
    sky = init * model[:, np.newaxis]
    sky[~selg] = 0.0
    res = 0.0 * sky
    skysub = rect - sky
    for j in np.arange(rect.shape[1]):
        E = error[:, j] * 1.0
        E[E == 0.0] = 1e9
        W = 1.0 / E ** 2
        W = np.sqrt(np.diag(W))
        A = model[:, np.newaxis]
        B = skysub[:, j]
        Aw = np.dot(W, A)
        Bw = np.dot(B, W)
        solj = np.linalg.lstsq(Aw, Bw, rcond=None)[0][0]
        sg = np.sign(solj)
        V = [np.abs(solj), 0.5 * np.median(E)]
        mult = V[np.argmin(V)] * sg
        res[:, j] = mult
    return sky + res
