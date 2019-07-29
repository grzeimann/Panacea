#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:04:57 2019

@author: gregz
"""

import argparse as ap
import fnmatch
import glob
import numpy as np
import os.path as op
import sys
import tarfile
import warnings

from astropy.io import fits
from astropy.stats import biweight_location
from datetime import timedelta, datetime
from distutils.dir_util import mkpath
from input_utils import setup_logging
from scipy.interpolate import interp2d

warnings.filterwarnings("ignore")


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("date",
                    help='''Date for reduction''',
                    type=str, default=None)

parser.add_argument("-r", "--rootdir",
                    help='''Directory for raw data. 
                    (just before date folders)''',
                    type=str, default='/work/03946/hetdex/maverick')

parser.add_argument("-f", "--filename",
                    help='''Filename of the master fits file''',
                    type=str, default=None)

parser.add_argument("-k", "--kind", help='''twi or flt''',  type=str,
                    default='twi')

parser.add_argument("-o", "--outputdir",
                    help='''Output Directory''',
                    type=str, default='powerlaw')

args = parser.parse_args(args=None)

log = setup_logging()
args.log = log

if args.kind not in ['zro', 'drk', 'sci', 'twi', 'flt']:
    args.log.error('"--kind" argument did not match "zro" or "drk"')
    sys.exit(1)

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
    plaw = I(xind[0, :], yind[:, 0]).swapaxes(0, 1)
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


def get_twi_spectra(array_flt, array_trace):
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
    twi_spectrum = np.zeros((array_trace.shape[0], array_trace.shape[1]))
    N = array_flt.shape[0]
    x = np.arange(array_flt.shape[1])
    for fiber in np.arange(array_trace.shape[0]):
        if array_trace[fiber].min() < 1.:
            continue
        if np.ceil(array_trace[fiber]).max() >= (N-1):
            continue
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        twi_spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
    twi_spectrum[~np.isfinite(twi_spectrum)] = 0.0
    return twi_spectrum

def splitall(path):
    ''' 
    Split path into list 
    
    Parameters
    ----------
    path : str
        file path
    
    Returns
    -------
    allparts : str
        list of constituent parts of the path
    '''
    allparts = []
    while 1:
        parts = op.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def get_ifuinfo_from_header(tfile, fn):
    ''' 
    Grab RA, Dec, and PA from header
    The fits header may be in a tar file, "tfile"
    
    Parameters
    ----------
    tfile : str
        Tarfile name
    fn : str
        fits filename containing the header in the primary HDU
    
    Returns
    -------
    specid : str
    ifuid : str
    ifuslot : str
    '''
    if tfile is not None:
        t = tarfile.open(tfile,'r')
        a = fits.open(t.extractfile(fn))
    else:
        a = fits.open(fn)
    specid = '%03d' % int(a[0].header['SPECID'])
    ifuid = '%03d' % int(a[0].header['IFUID'])
    ifuslot = '%03d' % int(a[0].header['IFUSLOT'])
    amp = (a[0].header['CCDPOS'].replace(' ', '') +
           a[0].header['CCDHALF'].replace(' ', ''))
    return specid, ifuid, ifuslot, a[0].header, amp


def build_path(rootdir, date, obs, ifuslot, amp, base='sci', exp='exp*',
               instrument='virus'):
    '''
    Build path for a given ifuslot, amplifier, observation, date, and rootdir
    
    Parameters
    ----------
    rootdir : str
        root directory where raw data is stored
    date : str
        date for path
    obs : integer
        observation for path
    ifuslot : str
        IFU slot for path
    amp : str
        Amplifier for path
    base : str
        Type of exposure
    exp : str
        Exposure number for path
    instrument : str
        Instrument for path
    
    Returns
    -------
    path : str
        File path constructed for the given date, obs, ect.
    '''
    if obs != '*':
        obs = '%07d' % obs
    path = op.join(rootdir, date, instrument, '%s%s' % (instrument, obs),
                   exp, instrument, '2*%s%s*%s.fits' % (ifuslot, amp, base))
    return path
    

def get_ifuslots(filenames):
    '''
    Get ifuslots for the rootdir, date, and observation in the global variable
    args.
    
    Parameters
    ----------
    filenames : list
        list of filenames
    
    Returns
    -------
    ifuslots : 1d numpy array [type=int]
        ifuslots for the input date, observation
    '''
    ifuslots = [int(op.basename(fn).split('_')[1][:3]) for fn in filenames]
    return np.unique(ifuslots)

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


def base_reduction(filename, get_header=False, tfile=None):
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
    if tfile is not None:
        t = tarfile.open(tfile,'r')
        a = fits.open(t.extractfile(filename))
    else:
        a = fits.open(filename)

    image = np.array(a[0].data, dtype=float)
    
    # Overscan subtraction
    overscan_length = 32 * (image.shape[1] / 1064)
    O = biweight_location(image[:, -(overscan_length-2):])
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


def get_mastertwi(files, masterbias, twitarfile):
    '''
    Make a master flat image from the twilight frames
    
    Parameters
    ----------
    files : list
        list of fits file names
    masterbias : 2d numpy array
        masterbias in e-
    twitarfile : str or None
        Name of the tar file if the fits file is tarred
    
    Returns
    -------
    mastertwi : 2d numy array
        median stacked twilight frame
    '''
    listtwi = []
    for filename in files:
        a, e = base_reduction(filename, tfile=twitarfile)
        a[:] -= masterbias
        listtwi.append(a)
    twi_array = np.array(listtwi, dtype=float)
    norm = np.median(twi_array, axis=(1, 2))[:, np.newaxis, np.newaxis]
    return np.median(twi_array / norm, axis=0) * np.median(norm)

def get_twi_tarfile(pathname, date, kind='twi'):
    '''
    Go through many dates if necessary and find the tar file that contains
    twilight exposures.
    
    Parameters
    ----------
    pathname : str
        File path to look for the tar file containing twilight exposures
    date : str
        Date for initial search for twilight exposures
    
    Returns
    -------
    twitarfile : str
        Name of the tarred file containing twilight exposures
    '''
    datec = date
    datec_ = datetime(int(date[:4]), int(date[4:6]), int(date[6:]))
    daten_ = datec_ + timedelta(days=1)
    daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
    pathnamen = pathname[:]
    twitarfile = None
    while twitarfile is None:
        datec_ = datetime(int(daten[:4]), int(daten[4:6]), int(daten[6:]))
        daten_ = datec_ - timedelta(days=1)
        daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
        pathnamen = pathname.replace(datec, daten)
        tarfolders = glob.glob(pathnamen)
        for tarfolder in tarfolders:
            T = tarfile.open(tarfolder, 'r')
            flag = True
            while flag:
                a = T.next()
                try:
                    name = a.name
                except:
                    flag = False
                if name[-5:] == '.fits':
                    if name[-8:-5] == kind:
                        twitarfile = tarfolder
                    flag = False
            if twitarfile is not None:
                break
    return twitarfile


def roll_through_dates(pathname, date):
    '''
    Get appropriate cals on the closest night back 60 days in time
    
    Parameters
    ----------
    pathname : str
        File path to look for the tar file containing twilight exposures
    date : str
        Date for initial search for twilight exposures
    
    Returns
    -------
    twinames : list
        list of twilight fits file names
    '''
    datec = date
    daten = date
    cnt = 0
    pathnamen = pathname[:]
    while len(glob.glob(pathnamen)) == 0:
        datec_ = datetime(int(daten[:4]), int(daten[4:6]), int(daten[6:]))
        daten_ = datec_ - timedelta(days=1)
        daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
        pathnamen = pathname.replace(datec, daten)
        cnt += 1
        log.info(pathnamen)
        if cnt > 60:
            log.error('Could not find cal within 60 days.')
            break
    return glob.glob(pathnamen)


def get_script_path():
    '''
    Get script path, aka, where does Remedy live?
    '''
    return op.dirname(op.realpath(sys.argv[0]))


def get_twi_files(kind='twi'):
    '''
    Get the names of the fits files for the science frames of interest
    and the twilight frames of interest.  If the files are tarred,
    return the names of the tar files as well.
    
    Parameters
    ----------
    args : object [GLOBAL]
        Object containing the command line inputs from the python call
    
    Returns
    -------
    twinames : list
        names of the twilight fits files
    twitarfile : str or None
        name of the tar file containing the fits files if there is one
    '''
    file_glob = build_path(args.rootdir, args.date, '*',
                           '*', 'LL')
    path = splitall(file_glob)
    tarname = op.join(*path[:-3]) + ".tar"
    if len(glob.glob(tarname)):
        twitarfile = get_twi_tarfile(tarname, args.date, kind=kind)
        with tarfile.open(twitarfile) as tf:
            twinames = [fn for fn in sorted(tf.getnames())
                        if fn[-5:] == '.fits']
    else:
        pathname = build_path(args.rootdir, args.date, '*', '*', '*',
                             base=kind)
        twinames = sorted(roll_through_dates(pathname, args.date))
        twitarfile = None
    return twinames, twitarfile

def main_all():
    twinames, twitarfile = get_twi_files(kind=args.kind)
    ifuslots = get_ifuslots(twinames)
    outdir = op.join(args.outputdir, args.date)
    mkpath(outdir)
    for ifuslot in ifuslots:
        for amp in ['LL', 'LU', 'RL', 'RU']:
            ifuSLOT = '%03d' % int(ifuslot)
            fnames_glob = '*/2*%s%s*%s.fits' % (ifuSLOT, amp, args.kind)
            twibase = fnmatch.filter(twinames, fnames_glob)
            log.info('Making powerlaw for %s%s' % (ifuslot, amp))
            try:
                masterbias = 0.0
                masterflt = get_mastertwi(twibase, masterbias, twitarfile)
                sinfo = get_ifuinfo_from_header(twitarfile, twibase[0])
                specid, ifuid, ifuSLOT, header, amP = sinfo
                trace, ref = get_trace(masterflt, specid, ifuSLOT, ifuid, amp,
                                       args.date)
                twi = get_twi_spectra(masterflt, trace)
                medtwi = np.median(twi, axis=0)
                plaw = get_powerlaw(masterflt, trace, twi, amp)
                plaw /= masterflt.sum()
                plaw *= masterflt.shape[0] * masterflt.shape[1]
                name = 'plaw_%s_%s_%s_%s.fits' % (specid, ifuSLOT, ifuid, amp)
                fits.PrimaryHDU(plaw, header=header).writeto(op.join(outdir,
                                name), overwrite=True)
            except:
                log.info('Failed to make powerlaw for %s%s' % (ifuslot, amp))


def main_single():
    outdir = op.join(args.outputdir, args.date)
    mkpath(outdir)
    masterflt = fits.open(args.filename)[0].data
    sinfo = get_ifuinfo_from_header(None, args.filename)
    specid, ifuid, ifuSLOT, header, amp = sinfo
    trace, ref = get_trace(masterflt, specid, ifuSLOT, ifuid, amp,
                           args.date)
    twi = get_twi_spectra(masterflt, trace)
    medtwi = np.median(twi, axis=0)
    plaw = get_powerlaw(masterflt, trace, twi, amp)
    plaw /= masterflt.sum()
    plaw *= masterflt.shape[0] * masterflt.shape[1]
    name = 'plaw_%s_%s_%s_%s.fits' % (specid, ifuSLOT, ifuid, amp)
    fits.PrimaryHDU(plaw, header=header).writeto(op.join(outdir,
                    name), overwrite=True)

if args.filename is not None:
    main_single()
else:
    main_all()
