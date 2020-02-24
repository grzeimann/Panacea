#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:36:42 2019

@author: gregz
"""

import glob
import numpy as np
import os.path as op
import tarfile

from astropy.io import fits
from math_utils import biweight
from datetime import datetime
from scipy.interpolate import interp2d, interp1d, LinearNDInterpolator

from astropy.modeling.models import Gaussian1D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clip, mad_std, sigma_clipped_stats
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.convolution import interpolate_replace_nans

from sklearn.cluster import AgglomerativeClustering


def orient_image(image, amp, ampname):
    '''
    Orient the images from blue to red (left to right)
    Fibers are oriented to match configuration files
    
    Parameters
    ----------
    image : 2d numpy array
        fits image
    amp : str
        Amplifier for the fits image
    ampname : str
        Amplifier name is the location of the amplifier
    
    Returns
    -------
    image : 2d numpy array
        Oriented fits image correcting for what amplifier it comes from
        These flips are unique to the VIRUS/LRS2 amplifiers
    '''
    if amp == "LU":
        image[:] = image[::-1, ::-1]
    if amp == "RL":
        image[:] = image[::-1, ::-1]
    if ampname is not None:
        if ampname == 'LR' or ampname == 'UL':
            image[:] = image[:, ::-1]
    return image


def base_reduction(filename, get_header=False):
    '''
    Reduce filename from tarfile or fits file.
    
    Reduction steps include:
        1) Overscan subtraction
        2) Trim image
        3) Orientation
        4) Gain Multiplication
        5) Error propagation
    
    Parameters
    ----------
    filename : str
        Filename of the fits file
    get_header : boolean
        Flag to get and return the header
    tfile : str
        Tar filename if the fits file is in a tarred file
    
    Returns
    -------
    a : 2d numpy array
        Reduced fits image, see steps above
    e : 2d numpy array
        Associated error frame
    '''
    # Load fits file
    tarbase = op.dirname(op.dirname(op.dirname(filename))) + '.tar'
    if op.exists(tarbase):
        T = tarfile.open(tarbase, 'r')
        s = '/'.join(filename.split('/')[-4:])
        a = fits.open(T.extractfile(s))
    else:
        a = fits.open(filename)

    image = np.array(a[0].data, dtype=float)
    
    # Overscan subtraction
    overscan_length = int(32 * (image.shape[1] / 1064))
    O = biweight(image[:, -(overscan_length-2):])
    image[:] = image - O
    
    # Trim image
    image = image[:, :-overscan_length]
    
    # Gain multiplication (catch negative cases)
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
    header = a[0].header
    
    # Orient image
    a = orient_image(image, amp, ampname) * gain
    
    # Calculate error frame
    E = np.sqrt(rdnoise**2 + np.where(a > 0., a, 0.))
    if get_header:
        return a, E, header
    return a, E

def power_law(x, c1, c2=.5, c3=.15, c4=1., sig=2.5):
    '''
    Power law for scattered light from mirror imperfections
    
    Parameters
    ----------
    x : float or 1d numpy array
        Distance from fiber in pixels
    c1 : float
        Normalization of powerlaw
    
    Returns
    -------
    plaw : float or 1d numpy array
        see function form below
    '''
    return c1 / (c2 + c3 * np.power(abs(x / sig), c4))


def get_powerlaw_ydir(trace, spec, amp, col):
    '''
    Get powerlaw in ydir for a given column
    
    Parameters
    ----------
    trace : 2d numpy array
        y position as function of x for each fiber
    spec : 2d numpy array
        fiber spectra
    amp : str
        amplifier
    col : int
        Column
    '''
    if amp in ['LL', 'RU']:
        ntrace = np.vstack([trace, 2064 - trace])
    else: 
        ntrace = np.vstack([-trace, trace])
    nspec = np.vstack([spec, spec])
    YM, XM = np.indices(ntrace.shape)
    yz = np.linspace(0, 1031, 25)
    plaw = []
    for yi in yz:
        d = np.sqrt((yi - ntrace[:, ::43])**2 + (col - XM[:, ::43])**2)
        plaw.append(np.nansum(nspec[:, ::43] *
                              power_law(d, 1.4e-5, c3=2., c4=1.0,  sig=1.5)))
    return yz, np.array(plaw)   


def get_powerlaw(image, trace, spec, amp):
    '''
    Solve for scatter light from powerlaw
    
    Parameters
    ----------
    image : 2d numpy array
        fits image
    trace : 2d numpy array
        y position as function of x for each fiber
    spec : 2d numpy array
        fiber spectra
    amp : str
        amplifier

    Returns
    -------
    plaw : 2d numpy array
        scatter light image from powerlaw
    '''
    fibgap = np.where(np.diff(trace[:, 400]) > 10.)[0]
    X = np.arange(image.shape[1])
    yind, xind = np.indices(image.shape)
    XV = np.array_split(X, 25)
    T = np.array_split(trace, 25, axis=1)
    XM, YM, ZM = ([], [], [])
    for xchunk, tchunk in zip(XV, T):
        avgy, avgz = ([], [])
        avgx = int(np.mean(xchunk))
        x, y = ([], [])
        dy = np.array(np.ceil(trace[0, xchunk])-7, dtype=int)
        for j, xc in enumerate(xchunk):
            d = np.arange(0, dy[j])
            if len(d):
                y.append(d)
                x.append([xc] * len(d))
        if len(y):
            y, x = [np.array(np.hstack(i), dtype=int) for i in [y, x]]
            avgy.append(np.mean(y))
            avgz.append(np.median(image[y, x]))
        for fib in fibgap:
            x, y = ([], [])
            dy = np.array(np.ceil(trace[fib, xchunk])+7, dtype=int)
            dy2 = np.array(np.ceil(trace[fib+1, xchunk])-7, dtype=int)
            for j, xc in enumerate(xchunk):
                d = np.arange(dy[j], dy2[j])
                if len(d):
                    y.append(d)
                    x.append([xc] * len(d))
            if len(y):
                y, x = [np.array(np.hstack(i), dtype=int) for i in [y, x]]
                avgy.append(np.mean(y))
                avgz.append(np.median(image[y, x]))
        x, y = ([], [])
        dy = np.array(np.ceil(trace[-1, xchunk])+7, dtype=int)
        for j, xc in enumerate(xchunk):
            d = np.arange(dy[j], image.shape[1])
            if len(d):
                y.append(d)
                x.append([xc] * len(d))
        if len(y):
            y, x = [np.array(np.hstack(i), dtype=int) for i in [y, x]]
            avgy.append(np.mean(y))
            avgz.append(np.median(image[y, x]))
        yz, plaw_col = get_powerlaw_ydir(trace, spec, amp, avgx)
        norm = np.nanmedian(np.array(avgz) / np.interp(avgy, yz, plaw_col))
        XM.append([avgx] * len(yz))
        YM.append(yz)
        ZM.append(plaw_col * norm)
    XM, YM = (np.hstack(XM), np.hstack(YM))
    xi, yi = (np.unique(XM), np.unique(YM))
    I = interp2d(xi, yi, np.hstack(ZM).reshape(len(yi), len(xi)), kind='cubic',
                 bounds_error=False)
    plaw = I(xind[0, :], yind[:, 0]).swapaxes(0, 1)[:, ::-1]
    return plaw


def get_trace_reference(specid, ifuslot, ifuid, amp, obsdate,
                       virusconfig='/work/03946/hetdex/maverick/virus_config'):
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


def get_trace(twilight, specid, ifuslot, ifuid, amp, obsdate, tr_folder):
    ref = get_trace_reference(specid, ifuslot, ifuid, amp, obsdate,
                              virusconfig=tr_folder)
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
        if len(trace) == len(ref):
            T = trace
        Trace[:, k] = T
        k += 1
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] > 0.
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)
    return trace, ref


def get_spectra(array_flt, array_trace, array_mod=None, npix=5):
    '''
    Extract spectra by dividing the flat field and averaging the central
    two pixels
    
    Parameters
    ----------
    array_flt : 2d numpy array
        twilight image
    array_trace : 2d numpy array
        trace for each fiber
    wave : 2d numpy array
        wavelength for each fiber
    def_wave : 1d numpy array [GLOBAL]
        rectified wavelength
    
    Returns
    -------
    twi_spectrum : 2d numpy array
        rectified twilight spectrum for each fiber  
    '''
    if array_mod is not None:
        pois_var = array_flt
        pois_var[pois_var<0.] = 0.
        error_array = np.sqrt(pois_var + 3.5**2)
        chi2 = np.zeros((array_trace.shape[0], array_trace.shape[1]))
        error = np.zeros((array_trace.shape[0], array_trace.shape[1]))
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))
    N = array_flt.shape[0]
    x = np.arange(array_flt.shape[1])
    LB = int((npix+1)/2)
    HB = -LB + npix + 1
    for fiber in np.arange(array_trace.shape[0]):
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N-LB):
            continue
        indv = np.round(array_trace[fiber]).astype(int)
        if array_mod is not None:
            chi2_a = np.zeros((3, array_trace.shape[1], HB+LB))
        for j in np.arange(-LB, HB):
            if j == -LB:
                w = indv + j + 1 - (array_trace[fiber] - npix/2.)
            elif j == HB-1:
                w = (npix/2. + array_trace[fiber]) - (indv + j) 
            else:
                w = 1.
            if array_mod is not None:
                chi2_a[0, :, j+LB] = array_flt[indv+j, x] * w
                chi2_a[1, :, j+LB] = array_mod[indv+j, x] * w
                chi2_a[2, :, j+LB] = error_array[indv+j, x] * w
                error[fiber] += error_array[indv+j, x]**2 * w
            spec[fiber] += array_flt[indv+j, x] * w
            
        if array_mod is not None:
            norm = chi2_a[0].sum(axis=1) / chi2_a[1].sum(axis=1)
            num = (chi2_a[0] - chi2_a[1] * norm[:, np.newaxis])**2
            denom = (chi2_a[2] + 0.01*chi2_a[0].sum(axis=1)[:, np.newaxis])**2
            chi2[fiber] = 1. / (1. + 5.) * np.sum(num / denom, axis=1)
    if array_mod is not None:
        return spec, chi2, np.sqrt(error)
    return spec


def get_ifucenfile(folder, ifuid, amp):
    if ifuid == '004':
        ifucen = np.loadtxt(op.join(folder, 'IFUcen_files',
                            'IFUcen_HETDEX_reverse_R.txt'),
                            usecols=[0,1,2,4], skiprows=30)
        ifucen[224:,:] = ifucen[-1:223:-1,:]
    else:
        ifucen = np.loadtxt(op.join(folder, 'IFUcen_files',
                                    'IFUcen_HETDEX.txt'),
                            usecols=[0,1,2,4], skiprows=30)
        if ifuid in ['003','005','008']:
            ifucen[224:,:] = ifucen[-1:223:-1,:]
        if ifuid == '007':
            a = ifucen[37,1:3] * 1.
            b = ifucen[38,1:3] * 1.
            ifucen[37,1:3] = b
            ifucen[38,1:3] = a
        if ifuid == '025':
            a = ifucen[208,1:3] * 1.
            b = ifucen[213,1:3] * 1.
            ifucen[208,1:3] = b
            ifucen[213,1:3] = a
        if ifuid == '030':
            a = ifucen[445,1:3] * 1.
            b = ifucen[446,1:3] * 1.
            ifucen[445,1:3] = b
            ifucen[446,1:3] = a
        if ifuid == '038':
            a = ifucen[302,1:3] * 1.
            b = ifucen[303,1:3] * 1.
            ifucen[302,1:3] = b
            ifucen[303,1:3] = a
        if ifuid == '041':
            a = ifucen[251,1:3] * 1.
            b = ifucen[252,1:3] * 1.
            ifucen[251,1:3] = b
            ifucen[252,1:3] = a

    if amp=="LL":
        return ifucen[112:224,1:3][::-1,:]
    if amp=="LU":
        return ifucen[:112,1:3][::-1,:]
    if amp=="RL":
        return ifucen[224:336,1:3][::-1,:]               
    if amp=="RU":
        return ifucen[336:,1:3][::-1,:] 
    
def safe_sigma_clip(y):
    try:
        mask = sigma_clip(y, masked=True, maxiters=None,
                          stdfunc=mad_std)
    except:
        mask = sigma_clip(y, iters=None, stdfunc=mad_std) 
    return mask

def identify_sky_pixels(sky, kernel=10.0):
    G = Gaussian1DKernel(kernel)
    cont = convolve(sky, G, boundary='extend')
    mask = safe_sigma_clip(sky - cont)
    for i in np.arange(5):
        nsky = sky * 1.
        mask.mask[1:] += mask.mask[:-1]
        mask.mask[:-1] += mask.mask[1:]
        nsky[mask.mask] = np.nan
        cont = convolve(nsky, G, boundary='extend')
        while np.isnan(cont).sum():
            cont = interpolate_replace_nans(cont, G)
        mask = safe_sigma_clip(sky - cont)
    return mask.mask, cont

def find_peaks(y, thresh=10.):
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
    loc = loc[peaks > thresh]+1
    peak_loc = get_peaks(y, loc)
    peaks = y[np.round(peak_loc).astype(int)]
    return peak_loc, peaks

def get_wave_single(y, T_array, order=3, thresh=5., dthresh=25.):
    mask, cont = identify_sky_pixels(y, kernel=1.5)
    bl, bm = biweight(y - cont, calc_std=True)
    thr = bm * thresh
    loc, peaks = find_peaks(y - cont, thresh=thr)
    x = np.arange(len(y))
    dw = T_array[-1] - T_array[0]
    dx = loc[-1] - loc[0]
    p = np.array([-9.83655146e-09,  1.18442454e-04, -4.58877122e-01,
                  5.75772866e+02])
    wave = (x-loc[0]) * dw / dx + T_array[0]
    wave = wave + np.polyval(p, wave)
    guess_wave = np.interp(loc, x, wave)
    diff = np.abs(guess_wave[:, np.newaxis] - T_array)
    ind = np.argmin(diff, axis=0)
    sel = np.min(diff, axis=0) < 2.5
    P = np.polyfit(loc[ind][sel], T_array[sel], order)
    yv = np.polyval(P, loc[ind])
    res = np.std((T_array-yv)[sel])
    wave = np.polyval(P, x)
    return wave, res

def get_wave(spec, trace, T_array, res_lim=1., order=3):
    w = trace * 0.
    r = np.zeros((trace.shape[0],))
    xi = np.hstack([np.arange(2, trace.shape[0], 8), trace.shape[0]-3])
    failed = True
    thresh=20.
    while failed:
        for j in xi:
            try:
                S = np.nanmedian(spec[j-2:j+3], axis=0)
                wave, res = get_wave_single(S, T_array, order=order,
                                            thresh=thresh)
                w[j] = wave
                r[j] = res
            except:
                pass
        sel = (np.array(r) < res_lim) * (np.array(r) > 0.)
        if sel.sum() < 7:
            thresh -= 5.
        else:
            failed = False
        if thresh < 5.:
            return None
    wave = w * 0.
    xi = np.hstack([np.arange(0, trace.shape[1], 24), trace.shape[1]-1])
    for i in xi:
        x = trace[:, i]
        y = w[:, i]
        wave[:, i] = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    x = np.arange(wave.shape[1])
    for j in np.arange(wave.shape[0]):
        sel = wave[j] > 0.
        wave[j] = np.polyval(np.polyfit(x[sel], wave[j][sel], order), x)
    return wave

def get_pixelmask(dark):
    y1 = dark - np.median(dark, axis=0)[np.newaxis, :]
    y1 = y1 - np.median(y1, axis=1)[:, np.newaxis]
    m, m1, s = sigma_clipped_stats(y1)
    mask = np.abs(y1) > 5 * s
    yind, xind = np.where(y1 > 100)
    for yi, xi in zip(yind, xind):
        if yi + 2 < mask.shape[0]:
            if mask[yi+1, xi]:
                mask[:, xi] = True
    return np.array(mask, dtype=int)

def measure_fiber_profile(image, spec, trace, wave, def_wave, xmin=400,
                          xmax=600):
    ''' Measure the fiber profile for an image
    
    Parameters
    ----------
    image : 2d numpy array
        reduced image without scatter light or background counts
    spec : 2d numpy array
        fiber spectra from an aperture method for normalization
    trace : 2d numpy array
        fiber trace
    wave : 2d numpy array
        fiber wavelength
    def_wave : 1d numpy array
        standard wavelength of the spec array
        
    Returns
    -------
    init : interp1d model
        interpolation model to build 1d fiber profile for any give fiber
    '''
    yind, xind = np.indices(image.shape)
    
    profile = []
    for fibert, fiberw, fibers in zip(trace, wave, spec):
        ospec = interp1d(def_wave, fibers, kind='quadratic', fill_value=0.0,
                         bounds_error=False)(fiberw)
        indl = int(np.max([0, np.min(fibert)-10.]))
        indh = int(np.min([image.shape[0], np.max(fibert)+10.]))
        foff = yind[indl:indh, xmin:xmax] - fibert[np.newaxis, xmin:xmax]
        V = image[indl:indh, xmin:xmax] / ospec[np.newaxis, xmin:xmax]
        sel = np.abs(foff) <= 6.
        profile.append([foff[sel], V[sel]])
    imodel = []
    smodel = []
    for k in np.arange(8, 112, 16):
        allx = np.hstack([p[0] for j, p in enumerate(profile)
                          if np.abs(j-k)<=8])
        ally = np.hstack([p[1] for j, p in enumerate(profile)
                          if np.abs(j-k)<=8])
        inorder = np.argsort(allx)
        
        xbin = np.array([np.median(chunk) for chunk in np.array_split(allx[inorder], 25)])
        ybin = np.array([np.median(chunk) for chunk in np.array_split(ally[inorder], 25)])
        peak_loc, peaks = find_peaks(ybin, thresh=0.4)
        if len(peak_loc) != 1:
            continue
        peak_loc = np.interp(peak_loc, np.arange(len(xbin)), xbin)
        valley_loc, valley = find_peaks(1. - ybin, thresh=0.3)
        valley_loc = np.interp(valley_loc, np.arange(len(xbin)), xbin)
        if len(valley_loc) != 2:
            continue
        valley = (1. - valley) / 2.
        ud = xbin - peak_loc[0]
        yd = ybin
        sel = np.abs(ud) < 4.
        x = np.hstack([ud[sel], valley_loc[0]-peak_loc[0], valley_loc[1]-peak_loc[0], -5.0, 5.0])
        yp = np.hstack([yd[sel], valley[0], valley[1], 0., 0.])
        yp = yp[np.argsort(x)]
        x = np.sort(x)
        init = interp1d(x, yp, kind='quadratic', fill_value=0.0, bounds_error=False)
        shift = peak_loc[0]
        imodel.append(init)
        smodel.append(shift)
    x = np.linspace(-6., 6.)
    mod = []
    for model in imodel:
        try:
            mod.append(model(x))
        except:
            print('oops.')
    mod = np.median(mod, axis=0)
    init = interp1d(x, mod, kind='quadratic', fill_value=0.0,
                    bounds_error=False)

    return init, smodel


def build_model_image(init, image, trace, wave, spec, def_wave):
    yind, xind = np.indices(image.shape)
    model_image = image * 0.
    for fibert, fiberw, fibers in zip(trace, wave, spec):
        ospec = np.interp(fiberw, def_wave, fibers, left=0.0, right=0.0)
        indl = int(np.max([0, np.min(fibert)-10.]))
        indh = int(np.min([image.shape[0], np.max(fibert)+10.]))
        foff = yind[indl:indh, :] - fibert[np.newaxis, :]
        sel = np.abs(foff) <= 6.
        yi = yind[indl:indh, :][sel]
        xi = xind[indl:indh, :][sel]
        model_image[yi, xi] += (init(foff[sel]) *
                                (np.ones((indh-indl, 1)) * ospec[np.newaxis, :])[sel]) 
    return model_image


def get_continuum(spectra, nbins=25):
    '''
    Get continuum from sky-subtracted spectra

    Parameters
    ----------
    spectra : 2d numpy array
        sky-subtracted spectra for each fiber

    Returns
    -------
    cont : 2d numpy array
        continuum normalization
    '''
    a = np.array([biweight(f, axis=1) 
                  for f in np.array_split(spectra, nbins, axis=1)]).swapaxes(0, 1)
    x = np.array([np.mean(xi) 
                  for xi in np.array_split(np.arange(spectra.shape[1]), nbins)])
    cont = np.zeros(spectra.shape)
    X = np.arange(spectra.shape[1])
    for i, ai in enumerate(a):
        sel = np.isfinite(ai)
        if np.sum(sel)>nbins/2.:
            I = interp1d(x[sel], ai[sel], kind='quadratic',
                         fill_value=np.nan, bounds_error=False)
            cont[i] = I(X)
        else:
            cont[i] = 0.0
    return cont


def detect_sources(dx, dy, spec, err, mask, def_wave, psf, ran, scale,
                   spec_res=5.6, thresh=5.):
    '''
    Detection algorithm
    
    Parameters
    ----------
    dx : 1d numpy array
        delta_x or delta_ra in " for each fiber compared to a given 0, 0
    dy : 1d numpy array
        delta_y or delta_dec in " for each fiber compared to a given 0, 0
    spec : 2d numpy array
        spectrum (sky-subtracted) for each fiber
    err : 2d numpy array
        error spectrum for each fiber
    chi : 2d numpy array
        chi2 for each spectrum compared to a fiber profile when extracted
    ftf : 2d numpy array
        initial fiber to fiber for each fiber
    adj : 2d numpy array
        adjusted fiber to fiber for better sky subtraction
    mask : 2d numpy array
        masked spectrum values for bad fiber to fiber, bad pixels, or bad chi2
    seeing : float
        spatial seeing fwhm
    spec_res : float
        spectral resolution in pixels (2A) and refers to radius not fwhm
    
    Returns
    -------
    '''
    N1 = int((ran[1] - ran[0]) / scale) + 1
    N2 = int((ran[3] - ran[2]) / scale) + 1
    gridx, gridy = np.meshgrid(np.linspace(ran[0], ran[1], N1),
                               np.linspace(ran[2], ran[3], N2))
    T = np.array([psf[1].ravel(), psf[2].ravel()]).swapaxes(0, 1)
    I = LinearNDInterpolator(T, psf[0].ravel(), fill_value=0.0)
    psfpixscale = np.abs(psf[1][0, 1] - psf[1][0, 0])
    area = 0.75**2 * np.pi
    cube = np.zeros((gridx.shape[0], gridx.shape[1], len(def_wave)))
    errcube = np.zeros((gridx.shape[0], gridx.shape[1], len(def_wave)))
    origcube = np.zeros((gridx.shape[0], gridx.shape[1], len(def_wave)))
    origerrcube = np.zeros((gridx.shape[0], gridx.shape[1], len(def_wave)))
    mask = ~np.isfinite(spec)
    G = Gaussian1DKernel(spec_res/2.35/(def_wave[1]-def_wave[0]))
    cont = get_continuum(spec, nbins=25)
    S = spec - cont
    for i in np.arange(gridx.shape[0]):
        for j in np.arange(gridx.shape[1]):
            xg = gridx[i, j]
            yg = gridy[i, j]
            sel = np.where(np.sqrt((dx-xg)**2 + (dy-yg)**2)<=4.0)[0]
            weights = I(dx[sel]-xg, dy[sel]-yg) * area / psfpixscale**2
            imask = ~(mask[sel])
            X = S[sel]*1.
            X[mask[sel]] = 0.0
            Y = err[sel]*1.
            Y[mask[sel]] = 0.0
            origcube[i, j, :] = np.sum(weights[:, np.newaxis]*X*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0)
            origerrcube[i, j, :] = np.sqrt(np.sum(weights[:, np.newaxis]*Y**2*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0))
            w = np.sum(weights[:, np.newaxis] * imask, axis=0)
            cube[i, j, w<0.7] = np.nan
            errcube[i, j, w<0.7] = np.nan
            WS = convolve(origcube[i, j], G, preserve_nan=True, fill_value=0.0) / np.sum(G.array**2)
            WE = np.sqrt(convolve(origerrcube[i, j]**2, G, preserve_nan=True, fill_value=0.0) / np.sum(G.array**2))
            cube[i, j, :] = WS
            errcube[i, j, :] = WE
    Y = cube / errcube
    bl, bm = biweight(Y.ravel(), calc_std=True)
    test = Y > thresh * bm
    L = np.zeros((0, 3))
    if test.sum():
        ids = np.where(test)
        Z = Y[ids[0],ids[1],ids[2]]
        sid = np.argsort(Z)[::-1]
        ids_sorted = (ids[0][sid], ids[1][sid], ids[2][sid])
        z = np.array([gridx[ids_sorted[0], ids_sorted[1]]*3.,
                      gridy[ids_sorted[0], ids_sorted[1]]*3.,
                      def_wave[ids_sorted[2]]]).swapaxes(0, 1)
    
        clustering = AgglomerativeClustering(n_clusters=None, 
                                             compute_full_tree=True,
                                             distance_threshold=50,
                                             linkage='complete').fit(z)
        
        z = np.array([gridx[ids_sorted[0], ids_sorted[1]],
                      gridy[ids_sorted[0], ids_sorted[1]],
                      def_wave[ids_sorted[2]]]).swapaxes(0, 1)
        US = np.unique(clustering.labels_)
        L = np.zeros((len(US), 5))
        K = np.zeros((len(US), len(def_wave), 3))
        fitter = LevMarLSQFitter()
        for i, ui in enumerate(US):
            sel = clustering.labels_ == ui
            L[i, 0] = np.mean(z[sel, 0])
            L[i, 1] = np.mean(z[sel, 1])
            L[i, 2] = np.mean(z[sel, 2])
            dsel = np.sqrt((gridx - L[i, 0])**2 + (gridy - L[i, 1])**2) < 2.5
            wi = int(np.interp(L[i, 2], def_wave, np.arange(len(def_wave))))
            x = gridx[dsel]
            y = gridy[dsel]
            v = cube[:, :, wi][dsel]
            fsel = np.isfinite(v)
            xc = np.sum(x[fsel]*v[fsel]) / np.sum(v[fsel])
            yc = np.sum(y[fsel]*v[fsel]) / np.sum(v[fsel])
            sel = np.where(np.sqrt((dx-xc)**2 + (dy-yc)**2)<=4.0)[0]
            weights = I(dx[sel]-xc, dy[sel]-yc) * area / psfpixscale**2
            if weights.sum() < 0.7:
                L[i, :] = 0.0
                continue
            imask = ~(mask[sel])
            X = S[sel]*1.
            X[mask[sel]] = 0.0
            Y = err[sel]*1.
            Y[mask[sel]] = 0.0
            spatial_spec = np.sum(weights[:, np.newaxis]*X*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0)
            spatial_spec_err = np.sqrt(np.sum(weights[:, np.newaxis]*Y**2*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0))
            X = spec[sel]*1.
            X[mask[sel]] = 0.0
            Y = err[sel]*1.
            Y[mask[sel]] = 0.0
            spatial_spec_or = np.sum(weights[:, np.newaxis]*X*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0)
            spatial_spec_err_or = np.sqrt(np.sum(weights[:, np.newaxis]*Y**2*imask, axis=0) / np.sum(weights[:, np.newaxis]**2 * imask, axis=0))
            wsel = np.where(np.abs(def_wave - L[i, 2]) <= 8.)[0]
            if (~np.isfinite(spatial_spec[wsel])).sum() > 0.:
                L[i, :] = 0.0
                continue
            G = Gaussian1D(mean=L[i, 2], stddev=spec_res/2.35)
            G.stddev.bounds = (4./2.35, 8./2.35)
            G.mean.bounds = (L[i, 2] - 4., L[i, 2] + 4.)
            fit = fitter(G, def_wave[wsel], spatial_spec[wsel])
            wc = fit.mean.value
            csel = np.where(np.abs(def_wave - L[i, 2]) <= 10.)[0]
            chi2 = (1. / (len(csel) - 3.) *
                    np.sum((fit(def_wave[csel])-spatial_spec[csel])**2 /
                           spatial_spec_err[csel]**2))
            L[i, 0] = xc
            L[i, 1] = yc
            L[i, 2] = wc
            L[i, 3] = fit.stddev.value
            L[i, 4] = chi2
            K[i, :, 0] = spatial_spec_or
            K[i, :, 1] = spatial_spec_err_or
            K[i, :, 2] = fit(def_wave)
            
            
    return cube, errcube, origcube, origerrcube, L, K

def measure_contrast(image, spec, trace, xmin=0,
                          xmax=1032, low_per=5, high_per=95):
    ''' Measure the fiber profile for an image
    
    Parameters
    ----------
    image : 2d numpy array
        reduced image without scatter light or background counts
    spec : 2d numpy array
        fiber spectra from an aperture method for normalization
    trace : 2d numpy array
        fiber trace
    wave : 2d numpy array
        fiber wavelength
    def_wave : 1d numpy array
        standard wavelength of the spec array
        
    Returns
    -------
    init : interp1d model
        interpolation model to build 1d fiber profile for any give fiber
    '''
    yind, xind = np.indices(image.shape)
    
    contrast = []
    I = np.arange(2, 112, 4)
    for fibert, fibers in zip(trace[I], spec[I]):
        ospec = fibers
        indl = int(np.max([0, np.min(fibert)-10.]))
        indh = int(np.min([image.shape[0], np.max(fibert)+10.]))
        foff = yind[indl:indh, xmin:xmax] - fibert[np.newaxis, xmin:xmax]
        V = image[indl:indh, xmin:xmax] / ospec[np.newaxis, xmin:xmax]
        sel = np.abs(foff) <= 6.
        xI = xind[indl:indh, xmin:xmax]
        xI = xI[sel]
        ind_chunks = np.array_split(np.argsort(xI), 10)
        x = foff[sel]
        y = V[sel]
        for inds in ind_chunks:
            inorder = np.argsort(x)
            xbin = np.array([np.median(chunk) for chunk in np.array_split(x[inorder], 25)])
            ybin = np.array([np.median(chunk) for chunk in np.array_split(y[inorder], 25)])
            peak_loc, peaks = find_peaks(ybin, thresh=0.4)
            if len(peak_loc) != 1:
                continue
            peak_loc = np.interp(peak_loc, np.arange(len(xbin)), xbin)
            peak = np.interp(peak_loc, xbin, ybin)
            valley_loc, valley = find_peaks(1. - ybin, thresh=0.3)
            valley_loc = np.interp(valley_loc, np.arange(len(xbin)), xbin)
            if len(valley_loc) != 2:
                continue
            valleys = np.interp(valley_loc, xbin, ybin)
            valley = np.mean(valleys)
            contrast.append((peak - valley*2.) / (peak + valley*2.))
    if not len(contrast):
        return None, None, None
    return np.percentile(contrast, 50), np.percentile(contrast, low_per), np.percentile(contrast, high_per)
