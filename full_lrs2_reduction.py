# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
import numpy as np
import os.path as op
import glob
import sys
import warnings
import argparse as ap
from datetime import datetime

from astropy.io import fits
from astropy.table import Table
from utils import biweight_location
from scipy.signal import savgol_filter, medfilt2d
from scipy.ndimage.filters import percentile_filter
from scipy.interpolate import interp1d, interp2d, griddata
from input_utils import setup_logging
from astrometry import Astrometry
from astropy.stats import biweight_midvariance

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--date",
                    help='''Date for reduction''',
                    type=str, default='20181108')

parser.add_argument("-s", "--side",
                    help='''Blue or Red''',
                    type=str, default='red')

args = parser.parse_args(args=None)

blueinfo = [['BL', 'uv', '503_056_7001', [3640., 4640.], ['LL', 'LU'],
             [4350., 4375.], ['Hg_B', 'Cd-A_B', 'FeAr_R']],
            ['BR', 'orange', '503_056_7001',
             [4660., 6950.], ['RU', 'RL'], [6270., 6470.],
             ['Hg_B', 'Cd-A_B', 'FeAr_R']]]
redinfo = [['RL', 'red', '502_066_7002', [6450., 8400.], ['LL', 'LU'],
            [7225., 7425.], ['Hg_R', 'Cd-A_B', 'FeAr_R']],
           ['RR', 'farred', '502_066_7002',
            [8275., 10500.], ['RU', 'RL'], [9280., 9530.],
            ['Hg_R', 'Cd-A_B', 'FeAr_R']]]

fplane_file = '/work/03730/gregz/maverick/fplane.txt'
twi_date = args.date
sci_date = args.date

# FOR LRS2
instrument = 'lrs2'
if args.side.lower() == 'blue':
    info_side = blueinfo
else:
    info_side = redinfo

dither_pattern = np.zeros((50, 2))

log = setup_logging('panacea_quicklook')

baseraw = '/work/03946/hetdex/maverick'


sci_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp%s',
                   '%s', '2*_%sLL*sci.fits')
sciflt_path = op.join(baseraw, twi_date,  '%s', '%s%s', 'exp*',
                      '%s', '2*_%sLL_twi.fits')
cmp_path = op.join(baseraw, twi_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL_cmp.fits')
bias_path = op.join(baseraw, twi_date, '%s', '%s%s', 'exp*',
                    '%s', '2*_%sLL_zro.fits')


def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))


def get_ifucenfile(side, amp,
                   virusconfig='/work/03946/hetdex/maverick/virus_config',
                   skiprows=4):

    file_dict = {"uv": "LRS2_B_UV_mapping.txt",
                 "orange": "LRS2_B_OR_mapping.txt",
                 "red": "LRS2_R_NR_mapping.txt",
                 "farred": "LRS2_R_FR_mapping.txt"}

    ifucen = np.loadtxt(op.join(virusconfig, 'IFUcen_files',
                        file_dict[side]), usecols=[0, 1, 2], skiprows=skiprows)

    if amp == "LL":
        return ifucen[140:, 1:3][::-1, :]
    if amp == "LU":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RL":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RU":
        return ifucen[140:, 1:3][::-1, :]


def orient_image(image, amp, ampname):
    '''
    Orient the images from blue to red (left to right)
    Fibers are oriented to match configuration files
    '''
    if amp == "LU":
        image[:] = image[::-1, ::-1]
    if amp == "RL":
        image[:] = image[::-1, ::-1]
    if ampname is not None:
        if ampname == 'LR' or ampname == 'UL':
            image[:] = image[:, ::-1]
    return image


def make_avg_spec(wave, spec, binsize=35, per=50):
    ind = np.argsort(wave.ravel())
    T = 1
    for p in wave.shape:
        T *= p
    wchunks = np.array_split(wave.ravel()[ind],
                             T / binsize)
    schunks = np.array_split(spec.ravel()[ind],
                             T / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    nspec = np.array([np.percentile(chunk, per) for chunk in schunks])
    nwave, nind = np.unique(nwave, return_index=True)
    return nwave, nspec[nind]


def base_reduction(filename):
    a = fits.open(filename)
    image = np.array(a[0].data, dtype=float)
    # overscan sub
    overscan_length = 32 * (image.shape[1] / 1064)
    O = biweight_location(image[:, -(overscan_length-2):])
    image[:] = image - O
    # trim image
    image = image[:, :-overscan_length]
    gain = a[0].header['GAIN']
    gain = np.where(gain > 0., gain, 0.85)
    rdnoise = a[0].header['RDNOISE']
    rdnoise = np.where(rdnoise > 0., rdnoise, 3.)
    amp = (a[0].header['CCDPOS'].replace(' ', '') +
           a[0].header['CCDHALF'].replace(' ', ''))
    try:
        ampname = a[0].header['AMPNAME']
    except:
        ampname = None
    a = orient_image(image, amp, ampname) * gain
    E = np.sqrt(rdnoise**2 + np.where(a > 0., a, 0.))
    return a, E


def power_law(x, c1, c2=.5, c3=.15, c4=1., sig=2.5):
        return c1 / (c2 + c3 * np.power(abs(x / sig), c4))


def get_powerlaw(image, trace, spec):
    YM, XM = np.indices(image.shape)
    x = np.arange(image.shape[1])
    inds = []
    for j in np.arange(trace.shape[0]):
        inds.append(np.where(np.abs(YM - trace[j, np.newaxis, :]).ravel() <
                             7.)[0])
    inds = np.hstack(inds)
    inds = np.unique(inds)
    inds = np.setdiff1d(np.arange(image.shape[0] * image.shape[1]), inds)
    y, x = np.unravel_index(inds, image.shape)
    xlim = [0, image.shape[1]]
    xy = np.median(spec, axis=1)
    bottom = np.percentile(xy, 15)
    sel = np.where(xy > (3. * bottom))[0]
    xp = np.hstack([np.arange(xlim[0], xlim[1], 128), image.shape[1]])
    ylim = [0, image.shape[0]]
    yz = np.hstack([np.arange(ylim[0], ylim[1], 128), image.shape[0]])
    plaw, XX, YY = ([], [], [])
    YM, XM = np.indices(trace.shape)
    for xi in xp:
        for yi in yz:
            d = np.sqrt((yi - trace)**2 + (xi - XM)**2)
            plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2., c4=1.0,
                                                   sig=1.5)))
            XX.append(xi)
            YY.append(yi)
    for xi in xp:
        for s in sel:
            y0 = int(trace[s, xi])
            for i in np.arange(-6, 8, 2):
                d = np.sqrt(((y0+i) - trace)**2 + (xi - XM)**2)
                plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2.,
                                                       c4=1.0, sig=1.5)))
                YY.append(y0+i)
                XX.append(xi)
    plaw, XX, YY = [np.hstack(j) for j in [plaw, XX, YY]]
    grid_x, grid_y = np.meshgrid(np.arange(image.shape[1]),
                                 np.arange(image.shape[0]))
    C = griddata(np.array([XX, YY]).swapaxes(0, 1), plaw, (grid_x, grid_y),
                 method='cubic', fill_value=0.0).swapaxes(0, 1)
    norm = np.zeros((image.shape[1],))
    for b in np.arange(image.shape[1]):
        sel = (x == b)
        norm[b] = np.median(image[y[sel], x[sel]] / C[y[sel], x[sel]])
    return C * savgol_filter(norm, 141, 1)[np.newaxis, :]


def get_twiflat_field(flt_path, amp, array_wave, array_trace, common_wave,
                      masterbias, log):
    files = glob.glob(flt_path.replace('LL', amp))
    listflat = []
    array_flt, error = base_reduction(files[0])

    bigW = np.zeros(array_flt.shape)
    Y, X = np.indices(array_wave.shape)
    YY, XX = np.indices(array_flt.shape)
    for x, at, aw, xx, yy in zip(np.array_split(X, 2, axis=0),
                                 np.array_split(array_trace, 2, axis=0),
                                 np.array_split(array_wave, 2, axis=0),
                                 np.array_split(XX, 2, axis=0),
                                 np.array_split(YY, 2, axis=0)):
        for j in np.arange(at.shape[1]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p0 = np.polyfit(at[:, j], aw[:, j], 7)
            bigW[yy[:, j], j] = np.polyval(p0, yy[:, j])
    listspec = []
    for filename in files:
        log.info('Working on twiflat %s' % filename)
        array_flt, error = base_reduction(filename)
        array_flt[:] -= masterbias
        x = np.arange(array_wave.shape[1])
        spectrum = array_trace * 0.
        for fiber in np.arange(array_wave.shape[0]):
            indl = np.floor(array_trace[fiber]).astype(int)
            indh = np.ceil(array_trace[fiber]).astype(int)
            spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
        smooth = savgol_filter(spectrum, 315, 1, axis=1)
        avg = biweight_location(smooth, axis=(0,))
        norm = biweight_location(smooth / avg, axis=(1,))
        nw, ns = make_avg_spec(array_wave, spectrum / norm[:, np.newaxis],
                               binsize=41, per=95)
        I = interp1d(nw, ns, kind='linear', fill_value='extrapolate')
        ftf = spectrum * 0.
        for fiber in np.arange(array_wave.shape[0]):
            model = I(array_wave[fiber])
            ftf[fiber] = savgol_filter(spectrum[fiber] / model, 151, 1)
        nw1, ns1 = make_avg_spec(array_wave, spectrum / ftf, binsize=41,
                                 per=95)
        I = interp1d(nw1, ns1, kind='quadratic', fill_value='extrapolate')
        modelimage = I(bigW)
        flat = array_flt / modelimage
        listflat.append(flat)
        listspec.append(I(common_wave))
    flat = np.median(listflat, axis=(0,))
    flat[~np.isfinite(flat)] = 0.0
    flat[flat < 0.0] = 0.0
    return flat, bigW, np.nanmedian(listspec, axis=0)


def get_spectra(array_flt, array_trace):
    spectrum = array_trace * 0.
    x = np.arange(array_flt.shape[1])
    for fiber in np.arange(array_trace.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
    return spectrum


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


def find_cosmics(Y, E, thresh=8.):
    A = medfilt2d(Y, (5, 1))
    S = safe_division((Y - A), E)
    P = S - medfilt2d(S, (1, 15))
    x, y = np.where(P > thresh)
    xx, yy = ([], [])
    for i in np.arange(-0, 1):
        for j in np.arange(-0, 1):
            sel = ((x + i) >= 0) * ((x + i) < Y.shape[0])
            sel2 = ((y + j) >= 0) * ((y + j) < Y.shape[1])
            sel = sel * sel2
            xx.append((x + i)[sel])
            yy.append((y + j)[sel])
    xx = np.hstack(xx)
    yy = np.hstack(yy)
    inds = np.ravel_multi_index([xx, yy], Y.shape)
    inds = np.unique(inds)
    C = np.zeros(Y.shape, dtype=bool).ravel()
    C[inds] = True
    C = C.reshape(Y.shape)
    log.info('Number of pixels affected by cosmics: %i' % len(x))
    log.info('Fraction of pixels affected by cosmics: %0.5f' %
             (1.*len(inds)/Y.shape[0]/Y.shape[1]))
    return C


def weighted_extraction(image, error, flat, trace):
    E = safe_division(error, flat)
    E[E < 1e-8] = 1e9
    Y = safe_division(image, flat)
    cosmics = find_cosmics(Y, E)
    x = np.arange(trace.shape[1])
    spectrum = 0. * trace
    TT = np.zeros((trace.shape[0], 3, trace.shape[1], 4))
    for fiber in np.arange(trace.shape[0]):
        T = np.zeros((3, trace.shape[1], 4))
        indl = np.floor(trace[fiber]).astype(int)
        flag = False
        for ss, k in enumerate(np.arange(-1, 3)):
            try:
                T[0, :, ss] = Y[indl+k, x]
                T[1, :, ss] = 1. / E[indl+k, x]**2
                T[2, :, ss] = ~cosmics[indl+k, x]
            except:
                v = indl+k
                sel = np.where((v >= 0) * (v < Y.shape[0]))[0]
                T[0, sel, ss] = Y[v[sel], x[sel]]
                T[1, sel, ss] = 1. / E[v[sel], x[sel]]**2
                T[2, sel, ss] = ~cosmics[v[sel], x[sel]]
                flag = True
        if flag:
            if np.mean(indl) > (Y.shape[0]/2.):
                k = 2
            else:
                k = -1
            v = indl+k
            sel = np.where((v >= 0) * (v < len(x)))[0]
            a = np.sum(T[0, sel] * T[1, sel] * T[2, sel], axis=1)
            b = np.sum(T[1, sel] * T[2, sel], axis=1)
            spectrum[fiber, sel] = safe_division(a, b)
        else:
            a = np.sum(T[0] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            spectrum[fiber] = safe_division(a, b)
        TT[fiber] = T
    return spectrum


def get_trace_shift(sci_array, flat, array_trace, Yx):
    YM, XM = np.indices(flat.shape)
    inds = np.zeros((3, array_trace.shape[0], array_trace.shape[1]))
    XN = np.round(array_trace)
    inds[0] = XN - 1.
    inds[1] = XN + 0.
    inds[2] = XN + 1.
    inds = np.array(inds, dtype=int)
    Trace = array_trace * 0.
    FlatTrace = array_trace * 0.
    N = YM.max()
    x = np.arange(array_trace.shape[1])
    for i in np.arange(Trace.shape[0]):
        sel = YM[inds[0, i, :], x] >= 0.
        sel = sel * (YM[inds[2, i, :], x] < N)
        xmax = (YM[inds[1, i, sel], x[sel]] -
                (sci_array[inds[2, i, sel], x[sel]] -
                sci_array[inds[0, i, sel], x[sel]]) /
                (2. * (sci_array[inds[2, i, sel], x[sel]] -
                 2. * sci_array[inds[1, i, sel], x[sel]] +
                 sci_array[inds[0, i, sel], x[sel]])))
        Trace[i, sel] = xmax
        xmax = (YM[inds[1, i, sel], x[sel]] - (flat[inds[2, i, sel], x[sel]] -
                flat[inds[0, i, sel], x[sel]]) /
                (2. * (flat[inds[2, i, sel], x[sel]] - 2. *
                 flat[inds[1, i, sel], x[sel]] +
                 flat[inds[0, i, sel], x[sel]])))
        FlatTrace[i, sel] = xmax
    shifts = np.nanmedian(FlatTrace - Trace, axis=1)
    fits.PrimaryHDU(shifts).writeto('test_shifts.fits', overwrite=True)
    sys.exit(1)
    shifts = np.polyval(np.polyfit(np.nanmedian(FlatTrace, axis=1), shifts, 2),
                        Yx)
    return shifts


def extract_sci(sci_path, amps, flat, array_trace, array_wave, bigW,
                masterbias):
    files1 = sorted(glob.glob(sci_path.replace('LL', amps[0])))
    files2 = sorted(glob.glob(sci_path.replace('LL', amps[1])))

    array_list = []
    for filename1, filename2 in zip(files1, files2):
        log.info('Prepping sci %s' % filename1)
        array_flt1, e1 = base_reduction(filename1)
        array_flt2, e2 = base_reduction(filename2)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
    if len(array_list) > 1:
        sci_array = np.sum(array_list, axis=0)
    else:
        sci_array = np.squeeze(np.array(array_list))
    Xx = np.arange(flat.shape[1])
    Yx = np.arange(flat.shape[0])
    I = interp2d(Xx, Yx, flat, kind='cubic', bounds_error=False,
                 fill_value=0.0)
    shifts = get_trace_shift(sci_array, flat, array_trace, Yx)
    flat = I(Xx, Yx + shifts)
    log.info('Found shift for %s of %0.3f' % (files1[0], np.median(shifts)))
    array_list = []
    spec_list = []
    orig_list = []
    for filename1, filename2 in zip(files1, files2):
        log.info('Fiber extraction sci %s' % filename1)
        array_flt1, e1 = base_reduction(filename1)
        array_flt2, e2 = base_reduction(filename2)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_err = np.vstack([e1, e2])

        array_flt[:] -= masterbias
        array_list.append(array_flt)
        spectrum = weighted_extraction(array_flt, array_err, flat, array_trace)
        spectrum[~np.isfinite(spectrum)] = 0.0
        speclist = []
        for fiber in np.arange(array_wave.shape[0]):
            dlam = np.diff(array_wave[fiber])
            dlam = np.hstack([dlam[0], dlam])
            I = interp1d(array_wave[fiber], spectrum[fiber] / dlam,
                         kind='quadratic', fill_value='extrapolate')
            speclist.append(I(commonwave))
        spec_list.append(np.array(speclist))
        orig_list.append(spectrum)
    return np.array(array_list), np.array(spec_list), np.array(orig_list)


def get_masterbias(zro_path, amp):
    files = glob.glob(zro_path.replace('LL', amp))
    listzro = []
    for filename in files:
        a, error = base_reduction(filename)
        listzro.append(a)
    return np.median(listzro, axis=0)


def get_masterarc(arc_path, amp, arc_names, masterbias, specname):
    files = glob.glob(arc_path.replace('LL', amp))
    # if specname == 'farred':
    #     ofiles = glob.glob(arc_path.replace('LL', amp).replace('cmp', 'sci'))
    #     files = files + ofiles
    listarc, listarce = ([], [])
    for filename in files:
        f = fits.open(filename)
        # scicond = ('_056' in f[0].header['OBJECT']) and ('sci' in filename)
        if f[0].header['OBJECT'] in arc_names:
            a, e = base_reduction(filename)
            a[:] -= masterbias
            listarc.append(a)
            listarce.append(e)
    listarc, listarce = [np.array(x) for x in [listarc, listarce]]
    # total_error = np.sqrt(np.sum(listarce**2, axis=0))
    return np.sum(listarc, axis=0)


def get_mastertwi(twi_path, amp, masterbias):
    files = glob.glob(twi_path.replace('LL', amp))
    listtwi = []
    for filename in files:
        a, e = base_reduction(filename)
        a[:] -= masterbias
        listtwi.append(a)
    twi_array = np.array(listtwi, dtype=float)
    norm = np.median(twi_array, axis=(1, 2))[:, np.newaxis, np.newaxis]
    return np.median(twi_array / norm, axis=0)


def get_trace_reference(specid, ifuslot, ifuid, amp, obsdate,
                        virusconfig='/work/03946/hetdex/'
                        'maverick/virus_config'):
    files = glob.glob(op.join(virusconfig, 'Fiber_Locations', '*',
                              'fiber_loc_%s_%s_%s_%s.txt' %
                              (specid, ifuslot, ifuid, amp)))
    dates = [op.basename(op.dirname(fn)) for fn in files]
    obsdate = datetime(int(obsdate[:4]), int(obsdate[4:6]),
                       int(obsdate[6:]))
    timediff = np.zeros((len(dates),))
    for i, datei in enumerate(dates):
        d = datetime(int(datei[:4]), int(datei[4:6]),
                     int(datei[6:]))
        timediff[i] = np.abs((obsdate - d).days)
    ref_file = np.loadtxt(files[np.argmin(timediff)])
    return ref_file


def get_trace(twilight, specid, ifuslot, ifuid, amp, obsdate):
    ref = get_trace_reference(specid, ifuslot, ifuid, amp, obsdate)
    N1 = (ref[:, 1] == 0.).sum()
    good = np.where(ref[:, 1] == 0.)[0]

    def get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace
    image = twilight
    N = 40
    xchunks = np.array([np.mean(x) for x in
                        np.array_split(np.arange(image.shape[1]), N)])
    chunks = np.array_split(image, N, axis=1)
    flats = [np.median(chunk, axis=1) for chunk in chunks]
    Trace = np.zeros((len(ref), len(chunks)))
    k = 0
    for flat, x in zip(flats, xchunks):
        diff_array = flat[1:] - flat[:-1]
        loc = np.where((diff_array[:-1] > 0.) * (diff_array[1:] < 0.))[0]
        peaks = flat[loc+1]
        loc = loc[peaks > 0.1 * np.median(peaks)]+1
        trace = get_trace_chunk(flat, loc)
        T = np.zeros((len(ref)))
        if len(trace) == N1:
            T[good] = trace
            for missing in np.where(ref[:, 1] == 1)[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = (T[good[gind]] + ref[missing, 0] -
                              ref[good[gind], 0])
        Trace[:, k] = T
        k += 1
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] > 0.
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)
    return trace, ref


def find_peaks(y, thresh=8.):
    def get_peaks(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)
        Peaks = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Peaks
    diff_array = y[1:] - y[:-1]
    loc = np.where((diff_array[:-1] > 0.) * (diff_array[1:] < 0.))[0]
    peaks = y[loc+1]
    std = np.sqrt(biweight_midvariance(y))
    loc = loc[peaks > (thresh * std)]+1
    peak_loc = get_peaks(y, loc)
    peaks = y[np.round(peak_loc).astype(int)]
    return peak_loc, peaks/std


def robust_polyfit(x, y, order=3, niter=3):
    sel = y > 0.
    ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    for i in np.arange(niter):
        a = np.abs(y - ymod)
        mad = np.median(a)
        sel = a < 3. * mad
        if sel.sum() > (order+2):
            ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    return ymod


def count_matches(lines, loc, fib, cnt=5):
    found_lines = np.zeros((trace.shape[0], len(lines)))
    M = np.zeros((cnt, cnt))
    for k in np.arange(cnt):
        for j in np.arange(cnt):
            i1 = 0 + k
            i2 = -1 - j
            diff = [loc[fib][i1] - lines['col2'][0],
                    loc[fib][i2] - lines['col2'][-1]]
            m = (diff[1] - diff[0]) / (lines['col2'][-1] - lines['col2'][0])
            y = np.array(m * (lines['col2'] - lines['col2'][0]) +
                         diff[0] + lines['col2'])
            for i, line in enumerate(lines):
                col = y[i]
                v = np.abs(col - loc[fib])
                if np.min(v) < 3.:
                    found_lines[fib, i] = loc[fib][np.argmin(v)]
            M[k, j] = (found_lines[fib] > 0.).sum()
    return np.unravel_index(np.argmax(M), M.shape)


def get_wavelength_from_arc(image, trace, lines, side):
    if side == 'uv':
        thresh = 5.
    if side == 'orange':
        thresh = 50.
    if side == 'red':
        thresh = 50.
    if side == 'farred':
        thresh = 20.
    spectrum = get_spectra(image, trace)
    fib = np.argmax(np.median(spectrum, axis=1))
    cont = percentile_filter(spectrum, 15, (1, 101))
    spectrum -= cont
    x = np.arange(trace.shape[1])
    loc = []
    ph = []
    for i, spec in enumerate(spectrum):
        px, py = find_peaks(spec, thresh=thresh)
        loc.append(px)
        ph.append(py)
    ind1, ind2 = (0, 0)
    found_lines = np.zeros((trace.shape[0], len(lines)))
    diff = [loc[fib][ind1] - lines['col2'][0],
            loc[fib][-ind2-1] - lines['col2'][-1]]
    m = (diff[1] - diff[0]) / (lines['col2'][-1] - lines['col2'][0])
    y = np.array(m * (lines['col2'] - lines['col2'][0]) +
                 diff[0] + lines['col2'])
    # print(loc[fib], y, ph[fib])
    for i, line in enumerate(lines):
        col = y[i]
        v = np.abs(col - loc[fib])
        if np.min(v) < 5.:
            found_lines[fib, i] = loc[fib][np.argmin(v)]
    for i, line in enumerate(lines):
        if found_lines[fib, i] == 0.:
            continue
        for j in np.arange(0, fib)[::-1]:
            if len(loc[j]) < 1:
                continue
            k = j + 1
            v = found_lines[k, i]
            while (v == 0.) and (k < trace.shape[0]):
                k += 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.:
                found_lines[j, i] = loc[j][np.argmin(m)]
        for j in np.arange(fib+1, trace.shape[0]):
            if len(loc[j]) < 1:
                continue
            k = j - 1
            v = found_lines[k, i]
            while (v == 0.) and (k >= 0):
                k -= 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.:
                found_lines[j, i] = loc[j][np.argmin(m)]
    for i, line in enumerate(lines):
        if np.sum(found_lines[:, i]) < (0.5 * trace.shape[0]):
            found_lines[:, i] = 0.0
            continue
        ind = np.array(found_lines[:, i], dtype=int)
        xt = trace[np.arange(trace.shape[0]), ind]
        sel = found_lines[:, i] > 0.
        yt = robust_polyfit(xt, found_lines[:, i])
        found_lines[:, i] = yt
    wave = trace * 0.0
    res = np.zeros((trace.shape[0],))
    for j in np.arange(trace.shape[0]):
        sel = found_lines[j, :] > 0.0
        if sel.sum() < 3:
            continue
        wave[j] = np.polyval(np.polyfit(found_lines[j, sel],
                             lines['col1'][sel], 3), x)
        res[j] = np.std(np.interp(found_lines[j, sel], x, wave[j]) -
                        lines['col1'][sel])
    missing = np.where(np.all(wave == 0., axis=1))[0]
    if len(missing):
        good = np.where(~np.all(wave == 0., axis=1))[0]
        for j in np.arange(trace.shape[1]):
            x = trace[good, j]
            y = wave[good, j]
            wave[missing, j] = np.polyval(np.polyfit(x, y, 3),
                                          trace[missing, j])
    # print(res)
    log.info('Min, Max Wave: %0.2f, %0.2f' % (wave.min(), wave.max()))
    return wave


def get_objects(basefiles, attrs):
    s = []
    for fn in basefiles:
        F = fits.open(fn)
        s.append([])
        for att in attrs:
            s[-1].append(F[0].header[att])
    return s


def create_image_header(wave, xgrid, ygrid, zgrid, func=fits.ImageHDU):
    hdu = func(np.array(zgrid, dtype='float32'))
    hdu.header['CRVAL1'] = xgrid[0, 0]
    hdu.header['CRVAL2'] = ygrid[0, 0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CDELT1'] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header['CDELT2'] = ygrid[1, 0] - ygrid[0, 0]
    return hdu


def create_header_objection(wave, image, func=fits.ImageHDU):
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


def sky_subtraction(rect, xloc, yloc, seeing=1.5):
    D = np.sqrt((xloc[:, np.newaxis] - xloc)**2 + (yloc[:, np.newaxis] - yloc)**2)
    W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)
    N = W.sum(axis=0)

    def outlier(y, y1, oi):
        m = np.abs(y[oi] - y1[oi])
        o = (y - y1) > 3. * np.median(m)
        return o

    def fit_sky_col(x, y):
        #o = y == 0.
        #low = np.percentile(y[~o], 16)
        #mid = np.percentile(y[~o], 50)
        #high = np.percentile(y[~o], 84)
        smooth = (y[:, np.newaxis] * W).sum(axis=0) / N
        o = outlier(y, smooth, y>0.)
        for i in np.arange(3):
            smooth = (y[~o, np.newaxis] * W[~o]).sum(axis=0) / W[~o].sum(axis=0)
            o = outlier(y, smooth, ~o)
        return smooth

    x = np.arange(rect.shape[0])
    sky = rect * 0.
    for j in np.arange(rect.shape[1]):
        sky[:, j] = fit_sky_col(x, rect[:, j])
    return sky


def make_frame(xloc, yloc, data, wave, dw, Dx, Dy, wstart=5700.,
               wend=5800., scale=0.4, seeing_fac=1.3):
    seeing = seeing_fac * scale
    a, b = data.shape
    x = np.arange(xloc.min()-scale,
                  xloc.max()+1*scale, scale)
    y = np.arange(yloc.min()-scale,
                  yloc.max()+1*scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,)+xgrid.shape)
    area = 3. / 4. * np.sqrt(3.) * 0.59**2
    for k in np.arange(b):
        sel = np.isfinite(data[:, k])
        D = np.sqrt((xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid)**2 +
                    (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid)**2)
        W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)
        N = W.sum(axis=0)
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] *
                           W[sel]).sum(axis=0) / N * (scale**2 / area))
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')

    zimage = biweight_location(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


def write_cube(wave, xgrid, ygrid, zgrid, outname):
    hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
    hdu.header['CRVAL1'] = xgrid[0, 0]
    hdu.header['CRVAL2'] = ygrid[0, 0]
    hdu.header['CRVAL3'] = wave[0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CRPIX3'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CTYPE3'] = 'pixel'
    hdu.header['CDELT1'] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header['CDELT2'] = ygrid[1, 0] - ygrid[0, 0]
    hdu.header['CDELT3'] = wave[1] - wave[0]
    hdu.writeto(outname, overwrite=True)



# LRS2-R
fiberpos, fiberspec = ([], [])
log.info('Beginning the long haul.')
allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [],
                                                           [], [])

DIRNAME = get_script_path()

for info in redinfo:
    specinit, specname, multi, lims, amps, slims, arc_names = info
    arc_lines = Table.read(op.join(DIRNAME, 'lrs2_config/lines_%s.dat' %
                                   specname), format='ascii')
    commonwave = np.linspace(lims[0], lims[1], int(2064*1.5))
    specid, ifuslot, ifuid = multi.split('_')
    package = []
    for amp in amps:
        amppos = get_ifucenfile(specname, amp)
        ##############
        # MASTERBIAS #
        ##############
        log.info('Getting Masterbias for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        zro_path = bias_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        masterbias = get_masterbias(zro_path, amp)

        #####################
        # MASTERTWI [TRACE] #
        #####################
        log.info('Getting MasterTwi for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        twibase = sciflt_path % (instrument, instrument, '00000*', instrument,
                                 ifuslot)
        mastertwi = get_mastertwi(twibase, amp, masterbias)
        log.info('Getting Trace for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        trace, dead = get_trace(mastertwi, specid, ifuslot, ifuid, amp,
                                twi_date)
        fits.PrimaryHDU(trace).writeto('test_trace.fits', overwrite=True)

        ##########################
        # MASTERARC [WAVELENGTH] #
        ##########################
        log.info('Getting MasterArc for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        lamp_path = cmp_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        masterarc = get_masterarc(lamp_path, amp, arc_names, masterbias,
                                  specname)

        log.info('Getting Wavelength for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        wave = get_wavelength_from_arc(masterarc, trace, arc_lines, specname)
        fits.PrimaryHDU(wave).writeto('test_wave.fits', overwrite=True)

        #################################
        # TWILIGHT FLAT [FIBER PROFILE] #
        #################################
        log.info('Getting TwiFlat for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        twiflat, bigW, twispec = get_twiflat_field(twibase, amp, wave, trace,
                                                   commonwave, masterbias, log)
        package.append([wave, trace, twiflat, bigW, masterbias, amppos,
                        twispec, dead])
    # Normalize the two amps and correct the flat
    avg = package[0][-2] / 2. + package[1][-2] / 2.
    for i in np.arange(len(package)):
        norm = safe_division(package[i][-2], avg)
        I = interp1d(commonwave, norm, kind='quadratic',
                     fill_value='extrapolate')
        model = I(package[i][3])
        package[i][2][:] = package[i][2] * model
    calinfo = [np.vstack([package[0][i], package[1][i]])
               for i in np.arange(len(package[0]))]
    calinfo[1][package[0][1].shape[0]:, :] += package[0][2].shape[0]
    flatspec = get_spectra(calinfo[2], calinfo[1])
    calinfo.append(flatspec)
    # log.info('Getting Powerlaw of Flat Cal for %s' % specname)
    # plaw = get_powerlaw(calinfo[2], calinfo[1], flatspec)
    # fits.PrimaryHDU(plaw).writeto('test_plaw_%s.fits' % specname, overwrite=True)
    f = []
    for i, cal in enumerate(calinfo):
        if i == 0:
            func = fits.PrimaryHDU
        else:
            func = fits.ImageHDU
        f.append(func(cal))
    fits.HDUList(f).writeto('test_all_%s.fits' % specname, overwrite=True)
    #####################
    # SCIENCE REDUCTION #
    #####################
    basefiles = sorted(glob.glob(sci_path % (instrument, instrument, '0000*',
                                             '01', instrument, ifuslot)))
    all_sci_obs = [op.basename(op.dirname(op.dirname(op.dirname(fn))))[-7:]
                   for fn in basefiles]
    objects = get_objects(basefiles, ['OBJECT', 'EXPTIME'])
    for sci_obs, obj, bf in zip(all_sci_obs, objects, basefiles):
        log.info('Extracting %s from %s' % (obj[0], bf))
        scifiles = sci_path % (instrument, instrument, sci_obs, '*',
                               instrument, ifuslot)
        images, rect, spec = extract_sci(scifiles, amps, calinfo[2],
                                         calinfo[1], calinfo[0], calinfo[3],
                                         calinfo[4])
        cnt = 1
        wave_0 = np.mean(commonwave)
        darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
        T = Table.read(darfile, format='ascii.fixed_width_two_line')
        xoff = (np.interp(commonwave, T['wave'], T['x_0']) -
                np.interp(wave_0, T['wave'], T['x_0']))
        yoff = (np.interp(commonwave, T['wave'], T['y_0']) -
                np.interp(wave_0, T['wave'], T['y_0']))
        for im, r, s in zip(images, rect, spec):
            log.info('Subtracting sky %s, exp%02d' % (obj[0], cnt))
            r[calinfo[7][:, 1] == 1.] = 0.
            sky = sky_subtraction(r)
            skysub = r - sky
            outname = ('%s_%s_%s_%s_%s.fits' % ('multi', args.date, sci_obs,
                                                'exp%02d' % cnt, specname))
            cnt += 1
            X = np.array([T['wave'], T['x_0'], T['y_0']])
            for S, name in zip([sky, skysub], ['sky', 'skysub']):
                outname = ('%s_%s_%s_%s_%s_cube.fits' % (args.date, sci_obs,
                           'exp%02d' % cnt, specname, name))
                zcube, zimage, xgrid, ygrid = make_frame(calinfo[5][:, 0],
                                                         calinfo[5][:, 1], S,
                                                         commonwave,
                                                         T['wave'],
                                                         xoff, yoff,
                                                         wstart=wave_0-50.,
                                                         wend=wave_0+50.)
                write_cube(commonwave, xgrid, ygrid, zcube, outname)
            outname = ('%s_%s_%s_%s_%s.fits' % ('multi', args.date, sci_obs,
                                                'exp%02d' % cnt, specname))
            f1 = create_header_objection(commonwave, r, func=fits.PrimaryHDU)
            f2 = create_header_objection(commonwave, sky)
            f3 = create_header_objection(commonwave, skysub)
            f4 = create_image_header(commonwave, xgrid, ygrid, zimage)
            fits.HDUList([f1, f2, f3, f4, fits.ImageHDU(calinfo[5]),
                          fits.ImageHDU(commonwave),
                          fits.ImageHDU(X)]).writeto(outname, overwrite=True)

#        header = fits.open(glob.glob(sci_path % (instrument, instrument, sci_obs, '01',
#                                         instrument,
#                                         ifuslots[0]))[0])[0].header
#        PA = float(header['PARANGLE'])
#        RA = float(header['TRAJRA'])
#        DEC = float(header['TRAJDEC'])
#        log.info('Observation at %0.4f %0.4f, PA: %0.3f' % (RA, DEC, PA))
#        A = Astrometry(RA, DEC, PA, 0., 0., fplane_file=fplane_file)
#        for i in np.arange(nexp):
#            log.info('Getting RA, Dec for exposure, %i, ifuslot, %s, and amp,'
#                     ' %s' % (i+1, ifuslot, amp))
#            ra, dec = A.get_ifupos_ra_dec(ifuslot,
#                                          amppos[:, 0] + dither_pattern[i, 0],
#                                          amppos[:, 1] + dither_pattern[i, 1])
#            allra.append(ra)
#            alldec.append(dec)
#            allx.append(A.fplane.by_ifuslot(ifuslot).y + amppos[:, 0] +
#                        dither_pattern[i, 0])
#            ally.append(A.fplane.by_ifuslot(ifuslot).x + amppos[:, 1] +
#                        dither_pattern[i, 1])

#
#fitslist = [fits.PrimaryHDU(np.vstack(allspec)),
#            fits.ImageHDU(np.array(allflatspec)),
#            fits.ImageHDU(commonwave),
#            fits.ImageHDU(np.array(allra)),
#            fits.ImageHDU(np.array(alldec)),
#            fits.ImageHDU(np.array(allx)),
#            fits.ImageHDU(np.array(ally))]
#fits.HDUList(fitslist).writeto('test_big.fits', overwrite=True)
#flist1 = []
#alls = np.vstack(allsub)
#for j, resi in enumerate(alls):
#    if j == 0:
#        func = fits.PrimaryHDU
#    else:
#        func = fits.ImageHDU
#    flist1.append(func(resi))
#fits.HDUList(flist1).writeto('test_sub.fits', overwrite=True)
