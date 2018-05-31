#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""
import matplotlib
matplotlib.use('agg')
import argparse as ap
import numpy as np
import glob
import os.path as op
import os
import sys
from amplifier import Amplifier
from utils import biweight_location, biweight_midvariance
from CreateTexWriteup import CreateTex
from distutils.dir_util import mkpath
from astropy.io import fits
from operator import itemgetter
import logging
from scipy.signal import medfilt2d
import matplotlib.pyplot as plt
from fiber_utils import fit_continuum_sky, find_maxima
from utils import biweight_bin


matplotlib.rcParams['font.sans-serif'] = "Meiryo"
matplotlib.rcParams['font.family'] = "sans-serif"
# plt.style.use('seaborn-colorblind')

cmap = plt.get_cmap('Greys_r')

AMPS = ["LL", "LU", "RU", "RL"]


def write_fits(hdu, name):
    try:
        hdu.writeto(name, overwrite=True)
    except:
        hdu.writeto(name, clobber=True)


def setup_logging():
    '''Setup Logging for MCSED, which allows us to track status of calls and
    when errors/warnings occur.

    Returns
    -------
    log : class
        log.info() is for general print and log.error() is for raise cases
    '''
    log = logging.getLogger('characterize')
    if not len(log.handlers):
        # Set format for logger
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        fmt = logging.Formatter(fmt)
        # Set level of logging
        level = logging.INFO
        # Set handler for logging
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
        # Build log with name, mcsed
        log = logging.getLogger('characterize')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log


def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = '''Characterize a VIRUS calibration data set

    This script is used to characterized a set of calibration data,
    and can be run either on a dataset from the lab or a dataset from
    the mountain.

    (Note that the line breaks may require multiple copy and pastes)
    Example calls are as follows:

    python Panacea/characterize.py --rootdir '/data/characterization_lab'
    --output 'Characterized' -bd 20170909 -bo 3090010 -dd 20170909
    -do 3090012 -xd 20170909 -xo 3090005 -pd 20170909 -po 3090009
    -fd 20170908 -fo 3090001 --specid 309 --ifuslot 999 -q

    The description of each input parameter is listed below.'''

    parser = ap.ArgumentParser(description=description,
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument("--ifuslot", nargs='?', type=str,
                        help='''Single ifuslot value. [REQUIRED]
                        Ex: "075".''', default=None)
    parser.add_argument("-us", "--use_structure",
                        help='''Use well defined structure.''',
                        action="count", default=0)
    parser.add_argument("--date", nargs='?', type=str,
                        help='''If using "use_structure then [REQUIRED].
                        Ex: "20170912".''', default=None)
    parser.add_argument("--specid", nargs='?', type=str,
                        help='''Single specid value. [REQUIRED]
                        Ex: "304".''', default=None)
    parser.add_argument("--instr", nargs='?', type=str,
                        help='''Instrument to process.
                        Default: "camra"
                        Ex: "camra" for lab data,
                            "virus" for mountain.''', default="camra")
    parser.add_argument("--output", nargs='?', type=str,
                        help='''Output Directory
                        Default: \"characterized"''',
                        default="characterized")
    parser.add_argument("--rootdir", nargs='?', type=str,
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''',
                        default="/work/03946/hetdex/maverick")
    obstype = ['bia', 'drk', 'pxf', 'ptc', 'msk', 'ldl', 'arc']
    obsletter = ['b', 'd', 'x', 'p', 'm', 'l', 'a']
    obsname = ['Bias', 'Dark', 'Pixel Flat', 'Photon Transfer Curve',
               'Masked Fiber Flat', 'LDLS Fiber Flat', 'Arc Lamp']
    for t, l, n in zip(obstype, obsletter, obsname):
        parser.add_argument("-%sd" % l, "--%sdir_date" % t, nargs='?',
                            type=str, help=''' %s Directory Date.''' % n,
                            default=None)
        parser.add_argument("-%so" % l, "--%sdir_obsid" % t, nargs='?',
                            type=str, help=''' %s Directory Observation ID.'''
                            % n,
                            default=None)
        parser.add_argument("-%se" % l, "--%sdir_expnum" % t, nargs='?',
                            type=str, help=''' %s Directory Exposure Number.'''
                            % n,
                            default=None)
    parser.add_argument("-q", "--quick", help='''Quicker Version.''',
                        action="count", default=0)
    parser.add_argument("-dcb", "--dont_check_bias",
                        help='''Don't make masterbias.''',
                        action="count", default=0)
    parser.add_argument("-dcd", "--dont_check_dark",
                        help='''Don't make masterdark.''',
                        action="count", default=0)
    parser.add_argument("-dcr", "--dont_check_readnoise",
                        help='''Don't check the readnoise.''',
                        action="count", default=0)
    parser.add_argument("-dcg", "--dont_check_gain",
                        help='''Don't check the gain.''',
                        action="count", default=0)
    parser.add_argument("-dcp", "--dont_check_pixelflat",
                        help='''Don't make pixelflat.''',
                        action="count", default=0)
    parser.add_argument("-dcm", "--dont_check_mask",
                        help='''Don't check masked fiber flats''',
                        action="count", default=0)
    parser.add_argument("-dcl", "--dont_check_ldls",
                        help='''Don't check ldls fiber flats''',
                        action="count", default=0)
    args = parser.parse_args(args=argv)

    return args


def read_in_raw(args):
    log = setup_logging()
    # Check that the arguments are filled
    if args.ifuslot:
        args.ifuslot = "%03d" % int(args.ifuslot)
    else:
        msg = 'No IFUSLOT was provided, exiting now.'
        log.error(msg)
        sys.exit(1)

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations = []
    if args.use_structure:
        if args.date is None:
            msg = '"use_structure" is True but "--date" was not set.'
            msg += ' Exiting now.'
            log.error(msg)
            sys.exit(1)
        args.biadir_date = args.date
        args.biadir_obsid = '%03d%04d' % (int(args.specid), 10)
        args.drkdir_date = args.date
        args.drkdir_obsid = '%03d%04d' % (int(args.specid), 12)
        args.ptcdir_date = args.date
        args.ptcdir_obsid = '%03d%04d' % (int(args.specid), 9)
        args.pxfdir_date = args.date
        args.pxfdir_obsid = '%03d%04d' % (int(args.specid), 5)
        args.mskdir_date = args.date
        args.mskdir_obsid = '%03d%04d' % (int(args.specid), 3)
        args.ldldir_date = args.date
        args.ldldir_obsid = '%03d%04d' % (int(args.specid), 1)
        args.arcdir_date = args.date
        args.arcdir_obsid = '%03d%04d' % (int(args.specid), 2)
    if not args.dont_check_bias:
        observations.append('bia')
    if not args.dont_check_dark:
        observations.append('drk')
    if not args.dont_check_gain:
        observations.append('ptc')
    if not args.dont_check_pixelflat:
        observations.append('pxf')
    if not args.dont_check_mask:
        observations.append('msk')
    if not args.dont_check_ldls:
        observations.append('ldl')
        observations.append('arc')

    for obs in observations:
        amp_list = []
        for label in labels[:2]:
            getattr(args, obs+label)
            if getattr(args, obs+label) is None:
                msg = '%s%s was not provided.' % (obs, label)
                msg += ' Exiting now.'
                log.error(msg)
                sys.exit(1)
            else:
                setattr(args, obs+label,
                        getattr(args, obs+label).replace(" ", "").split(','))
        if getattr(args, obs+labels[2]) is not None:
            setattr(args, obs+labels[2],
                    getattr(args, obs+labels[2]).replace(" ", "").split(','))
        for date in getattr(args, obs+labels[0]):
            for obsid in getattr(args, obs+labels[1]):
                if getattr(args, obs+labels[2]) is not None:
                    for expnum in getattr(args, obs+labels[2]):
                        folder = op.join(date,
                                         args.instr,
                                         "{:s}{:07d}".format(args.instr,
                                                             int(obsid)),
                                         "exp{:02d}".format(int(expnum)),
                                         args.instr)
                        filepath = op.join(args.rootdir, folder,
                                           '*_%s*.fits' % args.ifuslot)
                        files = sorted(glob.glob(filepath))
                        if not len(files):
                            print('Found no files for path: %s' % filepath)
                        for fn in files:
                            amp = op.basename(fn).split('_')[1][-2:]
                            amp_list.append([fn, obs, amp])
                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr,
                                                         int(obsid)))
                    filepath = op.join(args.rootdir, folder, '*',
                                                     args.instr,
                                                     '*_%s*.fits'
                                                     % args.ifuslot)
                    files = sorted(glob.glob(filepath))
                    if not len(files):
                        print('Found no files for path: %s' % filepath)
                    for fn in files:
                        amp = op.basename(fn).split('_')[1][-2:]
                        amp_list.append([fn, obs, amp])

        setattr(args, obs + '_list', amp_list)

    return args


def make_plot(image_dict, outfile_name, vmin=-5, vmax=15):
    a, b = image_dict[AMPS[0]].shape
    fig = plt.figure(figsize=((1.*b/a)*4, 4))
    for i, amp in enumerate(AMPS):
        ax = plt.subplot(2, 2, i+1)
        ax.imshow(image_dict[amp], vmin=vmin, vmax=vmax, cmap=cmap,
                  origin='lower', interpolation='none')
        ax.text(b*.1, a*.7, amp, fontsize=24, color='r')
        ax.set_xticks([])
        ax.set_yticks([])
    plt.subplots_adjust(wspace=0.025, hspace=0.025)
    fig.savefig(outfile_name)


def make_ptc_plot(mn_dict, vr_dict, gain, rd, outfile_name, lowlim=100,
                  highlim=50000):
    fig = plt.figure(figsize=(6, 6))
    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,
                           figsize=(6, 6))
    xhl = np.log10(highlim)
    xll = np.log10(lowlim)
    yhl = np.log10(np.sqrt(highlim))
    yll = np.log10(np.sqrt(lowlim))
    cnt = 0
    x = np.logspace(xll, xhl)
    for i, row in enumerate(ax):
        for j, cell in enumerate(row):
            amp = AMPS[cnt]
            cell.plot(x, np.sqrt(1./gain[amp]*x), 'r', label='Shot')
            cell.plot(x, np.sqrt(rd[amp])*np.ones(x.shape), 'g',
                      label='Read Noise')
            cell.plot(x, np.sqrt(1./gain[amp]*x+rd[amp]), 'k',
                      label='Shot+Read')
            cell.plot(mn_dict[amp], np.sqrt(vr_dict[amp]), label='Measured')
            cell.text(10**(0.8*(xhl-xll)+xll), 10**(0.3*(yhl-yll)+yll), amp,
                      fontsize=24, color='r')
            cell.set_xlim([lowlim+0.5, highlim])
            cell.set_ylim([np.sqrt(lowlim)+1.5, np.sqrt(highlim)])
            cell.set_xscale('log')
            cell.set_yscale('log')
            if i == 0 and j == 0:
                cell.legend(loc='best', fancybox=True, framealpha=0.5)
            cnt += 1

    fig.text(0.5, 0.025, 'Signal', ha='center', fontsize=18)
    fig.text(0.025, 0.5, 'Noise', va='center', rotation='vertical',
             fontsize=18)
    plt.subplots_adjust(wspace=0.00, hspace=0.00)
    fig.savefig(outfile_name)


def check_bias(args, amp, folder, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    left_edge, right_edge, structure, overscan = [], [], [], []

    bia_list = []
    for itm in args.bia_list:
        if itm[2] == amp:
            bia_list.append(Amplifier(itm[0], '', name=itm[1]))
            bia_list[-1].subtract_overscan()
            bia_list[-1].trim_image()

    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i, v in enumerate(bia_list) if v.amp == amp]
    log = bia_list[sel[0]].log

    overscan_list = [[v.overscan_value for i, v in enumerate(bia_list)
                      if v.amp == amp]]
    overscan = biweight_location(overscan_list)
    log.info('Overscan value for %s: %0.3f' % (amp, overscan))
    # Loop through the bias list and measure the jump/structure
    big_array = np.array([v.image for v in itemgetter(*sel)(bia_list)])
    if args.quick:
        func = np.median
    else:
        func = biweight_location
    masterbias = func(big_array, axis=(0,))

    a, b = masterbias.shape
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))

    log.info('Writing masterbias_%s.fits' % (amp))
    write_fits(hdu,
               op.join(folder, 'masterbias_%s_%s.fits' % (args.specid, amp)))

    left_edge = func(masterbias[:, edge:edge+width])
    right_edge = func(masterbias[:, (b-width-edge):(b-edge)])
    structure = func(masterbias[:, edge:(b-edge)], axis=(0,))

    log.info('Left edge - Overscan, Right edge - Overscan: %0.3f, %0.3f'
             % (left_edge, right_edge))

    return left_edge, right_edge, structure, overscan, masterbias


def check_darks(args, amp, folder, masterbias, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    dark_counts = []
    drk_list = []
    for itm in args.drk_list:
        if itm[2] == amp:
            drk_list.append(Amplifier(itm[0], '', name=itm[1]))
            drk_list[-1].subtract_overscan()
            drk_list[-1].trim_image()

    # Select only the dark frames that match the input amp, e.g., "RU"
    sel = [i for i, v in enumerate(drk_list) if v.amp == amp]
    log = drk_list[sel[0]].log

    if len(sel) <= 2 or args.quick:
        func = np.median
    else:
        func = biweight_location
    log.info('Writing masterdark_%s.fits' % (amp))
    if len(sel) == 1:
        big_array = (v.image - masterbias)[np.newaxis, :, :]
    else:
        big_array = np.array([v.image - masterbias
                              for v in itemgetter(*sel)(drk_list)])
    masterdark = func(big_array, axis=(0,))
    a, b = masterdark.shape
    hdu = fits.PrimaryHDU(np.array(masterdark, dtype='float32'))
    write_fits(hdu,
               op.join(folder, 'masterdark_%s_%s.fits' % (args.specid, amp)))

    # Loop through the bias list and measure the jump/structure
    for s in sel:
        am = drk_list[s]
        a, b = am.image.shape
        dark_counts.append(func(am.image - masterbias) / am.exptime)
    s = biweight_location(dark_counts)
    log.info('Average Dark counts/s: %0.5f' % s)

    return s, masterdark


def measure_readnoise(args, amp):
    # Select only the bias frames that match the input amp, e.g., "RU"
    bia_list = []
    for itm in args.bia_list:
        if itm[2] == amp:
            bia_list.append(Amplifier(itm[0], '', name=itm[1]))
            bia_list[-1].subtract_overscan()
            bia_list[-1].trim_image()
    sel = [i for i, v in enumerate(bia_list) if v.amp == amp]
    log = bia_list[sel[0]].log
    # Make array of all bias images for given amp
    array_images = np.array([bia.image for bia in
                             itemgetter(*sel)(bia_list)])

    # Measure the biweight midvariance (sigma) for a given pixel and take
    # the biweight average over all sigma to reduce the noise in the first
    # measurement.
    if args.quick:
        func1 = np.median
        func2 = np.std
    else:
        func1 = biweight_location
        func2 = biweight_midvariance
    S = func1(func2(array_images, axis=(0,)))
    log.info("RDNOISE(ADU) for %s: %01.3f" % (amp, S))
    return S


def measure_gain(args, amp, rdnoise, flow=500, fhigh=35000, fnum=50):
    ptc_list = []
    for itm in args.ptc_list:
        if itm[2] == amp:
            ptc_list.append(Amplifier(itm[0], '', name=itm[1]))
            ptc_list[-1].subtract_overscan()
            ptc_list[-1].trim_image()
    sel = [i for i, v in enumerate(ptc_list) if v.amp == amp]
    log = ptc_list[sel[0]].log
    s_sel = list(np.array(sel)[
                 np.array([ptc_list[i].basename for i in sel]).argsort()])
    npairs = len(sel) / 2
    a, b = ptc_list[sel[0]].image.shape
    array_avg = np.zeros((npairs, a, b))
    array_diff = np.zeros((npairs, a, b))
    if args.quick:
        func1 = np.median
        func2 = np.std
    else:
        func1 = biweight_location
        func2 = biweight_midvariance
    for i in xrange(npairs):
        F1 = ptc_list[s_sel[2*i]].image
        F2 = ptc_list[s_sel[2*i+1]].image
        m1 = func1(F1)
        m2 = func1(F2)
        array_avg[i, :, :] = (F1 + F2) / 2.
        array_diff[i, :, :] = F1 * m2 / m1 - F2
    bins = np.logspace(np.log10(flow), np.log10(fhigh), fnum)
    gn = []
    array_avg = array_avg.ravel()
    array_diff = array_diff.ravel()
    mn_list = []
    vr_list = []
    for i in xrange(len(bins)-1):
        loc = np.where((array_avg > bins[i]) * (array_avg < bins[i+1]))[0]
        if len(loc) > 1e3:
            std = func2(array_diff[loc])
            vr = (std**2) / 2.
            vr_c = (std**2 - 2.*rdnoise**2) / 2.
            mn = func1(array_avg[loc])
            log.info("%s | Gain: %01.3f | RDNOISE (e-): %01.3f | <ADU>: %0.1f"
                     " | VAR: %0.1f | Pixels: %i"
                     % (amp, mn / vr_c, mn / vr_c * rdnoise, mn, vr, len(loc)))
            gn.append(mn / vr_c)
            mn_list.append(mn)
            vr_list.append(vr)
    sel = np.where((np.array(mn_list) > 1000.)*(np.array(mn_list) < 15000.))[0]
    if len(sel) > 2:
        s = func1(np.array(gn)[sel])
        log.info("Average Gain measurement for %s: %0.3f"
                 % (amp, s))
    else:
        log.warning("Not enough points for gain measurement, using -99.0")
        s = -99.
    return s, mn_list, vr_list, rdnoise**2


def make_pixelflats(args, amp, folder):
    pxf_list = []
    for itm in args.pxf_list:
        if itm[2] == amp:
            pxf_list.append(Amplifier(itm[0], '', name=itm[1]))
            pxf_list[-1].subtract_overscan()
            pxf_list[-1].trim_image()
    sel = [i for i, v in enumerate(pxf_list) if v.amp == amp]
    log = pxf_list[sel[0]].log

    a, b = pxf_list[sel[0]].image.shape
    masterflat = np.zeros((len(sel), a, b))

    for i, am in enumerate(itemgetter(*sel)(pxf_list)):
        masterflat[i, :, :] = am.image
    masterflat = np.median(masterflat, axis=(0,))
    smooth = medfilt2d(masterflat, (151, 1))
    masterflat = np.where(masterflat < 1e-8, 0.0, smooth / masterflat)
    smooth = medfilt2d(masterflat, (1, 151))
    pixflat = np.where(masterflat < 1e-8, 0.0, smooth / masterflat)

    hdu = fits.PrimaryHDU(np.array(pixflat, dtype='float32'))
    log.info('Writing pixelflat_%s.fits' % amp)
    write_fits(hdu, op.join(folder, 'pixelflat_%s.fits' % amp))

    return masterflat, pixflat


def power_law(x, c1, c2=.5, c3=.15, c4=1., sig=2.5):
        return c1 / (c2 + c3 * np.power(np.abs(x / sig), c4))


def make_master_image(args, amp_list, masterbias, masterdark, use_mean=False):
    ''' Make a master image from a selection in a list '''
    if len(amp_list) <= 2 or args.quick:
        if use_mean:
            func = np.mean
        else:
            func = np.median
    else:
        func = biweight_location
    big_array = np.array([v.image - masterbias - masterdark
                          for v in amp_list])
    master = func(big_array, axis=(0,))
    return master


def get_average_spec(fibers, nbins=1000):
    masterwave = []
    masterspec = []
    for fib, fiber in enumerate(fibers):
        masterwave.append(fiber.wavelength)
        masterspec.append(fiber.spectrum)
    masterwave = np.hstack(masterwave)
    masterspec = np.hstack(masterspec)
    nwave = np.linspace(masterwave.min(), masterwave.max(), nbins)
    return nwave, biweight_bin(nwave, masterwave, masterspec)


def check_ldls(args, amp, masterbias, masterdark, outname, folder, gain):
    ''' Works on contrast/fibermodel/wavelength/trace '''
    # Select only the bias frames that match the input amp, e.g., "RU"
    ldl_list = []
    for itm in args.ldl_list:
        if itm[2] == amp:
            ldl_list.append(Amplifier(itm[0], '', name=itm[1]))
            ldl_list[-1].subtract_overscan()
            ldl_list[-1].trim_image()
    sel = [i for i, v in enumerate(ldl_list) if v.amp == amp]
    log = ldl_list[sel[0]].log
    log.info('Writing masterflat_%s.fits' % (amp))
    masterflat = make_master_image(args, ldl_list, masterbias, masterdark)

    A = ldl_list[sel[0]]
    A.image = masterflat
    A.orient_image()
    hdu = fits.PrimaryHDU(np.array(A.image, dtype='float32'))
    write_fits(hdu, op.join(folder, 'masterflat_%s_%s.fits'
                                    % (args.specid, amp)))
    A.image_prepped = True
    A.use_trace_ref = False
    A.refit = True
    A.use_pixelflat = False
    A.gain = gain
    A.multiply_gain()
    A.check_fibermodel = True
    A.check_trace = False
    A.path = folder
    A.get_fibermodel()
    os.rename(op.join(folder, 'fibmodel_%s.png' % A.basename),
              op.join(folder, 'contrast_%s.png' % amp))
    A.fibers = get_wavelength_from_arc(args, amp, masterbias, masterdark,
                                       outname, folder, A.fibers)
    A.fiberextract()
    wave, avgspec = get_average_spec(A.fibers)
    waven = np.vstack([fiber.wavelength for fiber in A.fibers])
    specn = np.vstack([fiber.spectrum for fiber in A.fibers])
    waver, specr = rectify(waven, specn)
    hdu = fits.PrimaryHDU(np.array(specr, dtype='float32'))
    hdu.header['CRVAL1'] = waver[0]
    hdu.header['CDELT1'] = waver[1] - wave[0]
    write_fits(hdu, op.join(folder, 'Femasterflat_%s_%s.fits'
                                    % (args.specid, amp)))
    colors = plt.get_cmap('RdBu')(np.linspace(0., 1., len(A.fibers)))
    fig = plt.figure(figsize=(12, 8))
    for i, fiber in enumerate(A.fibers):
        plt.plot(fiber.wavelength, fiber.spectrum, color=colors[i],
                 alpha=0.3)
    plt.plot(wave, avgspec, color='magenta', lw=4, label='Average')
    plt.xlim([3480, 5530])
    plt.ylim([0., 300000.])
    plt.xlabel('Wavelength')
    plt.ylabel('e- per exposure')
    plt.legend()
    plt.savefig(op.join(folder, 'ldls_spectra_%s.png' % amp))
    plt.close(fig)
    return masterflat


def get_wavelength_from_arc(args, amp, masterbias, masterdark, outname, folder,
                            fibers):
    # Select only the bias frames that match the input amp, e.g., "RU"
    arc_list = []
    for itm in args.arc_list:
        if itm[2] == amp:
            arc_list.append(Amplifier(itm[0], '', name=itm[1]))
            arc_list[-1].subtract_overscan()
            arc_list[-1].trim_image()
    sel = [i for i, v in enumerate(arc_list) if v.amp == amp]
    log = arc_list[sel[0]].log
    log.info('Writing masterarc_%s.fits' % (amp))
    masterflat = make_master_image(args, arc_list, masterbias, masterdark,
                                   use_mean=True)

    A = arc_list[sel[0]]
    A.image = masterflat
    A.orient_image()
    A.image_prepped = True
    hdu = fits.PrimaryHDU(np.array(A.image, dtype='float32'))
    write_fits(hdu, op.join(folder, 'masterarc_%s_%s.fits'
                                    % (args.specid, amp)))
    A.fibers = list(fibers)
    A.fiberextract()
    wave_list = [[3652.1026, 78], [4046.5539, 277], [4077.8298, 293],
                 [4358.3253, 435], [4678.149, 596], [4799.912, 658],
                 [5085.822, 808], [5460.7366, 1005]]
    
    if len(A.fibers[0].spectrum) > 1032:
        thresh = 1e4
    else:
        thresh = 1e2
    for fiber in A.fibers:
        y = fiber.spectrum
        x = np.arange(len(y))
        d1 = np.diff(y)
        selu = np.where(d1 > thresh)[0]
        sell = np.where(d1 < -thresh)[0]
        ind = []
        for i in selu:
            cont = True
            for j in ind:
                if np.abs(j - i) < 5:
                    cont = False
            if cont:
                u = selu[np.where(np.abs(selu - i) < 10)[0]]
                l = sell[np.where(np.abs(sell - i) < 10)[0]]
                v = (u.sum() + l.sum()) / (len(u) + len(l))
                ind.append(v)
        fac = len(y) / 1032
        pr = np.array(ind) / fac
        d = []
        off = 0.0
        for wvi in wave_list:
            loc = np.argmin(np.abs(pr - wvi[1]))
            if np.abs(pr[loc] - wvi[1] - off) < 15*fac:
                off = pr[loc] - wvi[1]
                d.append([pr[loc]*fac, wvi[0]])
        d = np.array(d)
        p0 = np.polyfit(d[:, 0] / (len(y)*1.), d[:, 1], 3)
        fiber.wavelength = np.polyval(p0, x / (len(y)*1.))

    return A.fibers


def check_masked_fibers(args, amp, masterbias, masterdark, outname, folder):
    # Select only the bias frames that match the input amp, e.g., "RU"
    msk_list = []
    for itm in args.msk_list:
        if itm[2] == amp:
            msk_list.append(Amplifier(itm[0], '', name=itm[1]))
            msk_list[-1].subtract_overscan()
            msk_list[-1].trim_image()
    sel = [i for i, v in enumerate(msk_list) if v.amp == amp]
    log = msk_list[sel[0]].log
    log.info('Writing mastermaskflat_%s.fits' % (amp))
    mastermaskflat = make_master_image(args, msk_list, masterbias, masterdark)

    A = msk_list[sel[0]]
    A.image = mastermaskflat
    A.orient_image()
    hdu = fits.PrimaryHDU(np.array(A.image, dtype='float32'))
    write_fits(hdu, op.join(folder, 'mastermaskedflat_%s_%s.fits'
                                    % (args.specid, amp)))
    A.image_prepped = True
    A.use_trace_ref = False
    A.refit = True
    A.use_pixelflat = False
    A.trace_y_window = 50.
    A.trace_repeat_length = 40
    A.gain = 1.
    A.check_trace = False
    A.get_trace()

    n, d = A.image.shape
    col = np.arange(d)
    nwave = 3
    fsize = 15
    radius = 5.
    fibs = [2, 5]
    cols = np.arange(d)
    f, ax = plt.subplots(len(fibs), nwave, sharey=True, sharex=True,
                         figsize=(nwave*4, len(fibs)*4))
    stot = 0
    for fiber in A.fibers:
        llim = np.array(np.max([np.zeros((d,)), fiber.trace-radius], axis=0),
                        dtype=int)
        ulim = np.array(np.min([np.ones((d,))*n, fiber.trace+radius+1],
                               axis=0),
                        dtype=int)
        for ll, ul, c in zip(llim, ulim, cols):
            stot += A.image[ll:ul, c].sum()
    sbig = A.image.sum()
    splaw = 100. * (sbig - stot) / (sbig * 1.)
    f.suptitle('The percentage of flux in the powerlaw is: %0.3f%%' % splaw)
    for fi, fib in enumerate(fibs):
        fiber = A.fibers[fib]
        llim = np.array(np.max([np.zeros((d,)), fiber.trace-fsize], axis=0),
                        dtype=int)
        ulim = np.array(np.min([np.ones((d,))*n, fiber.trace+fsize+1], axis=0),
                        dtype=int)
        nstep = int((n-100.) / nwave)
        for i in np.arange(nwave):
            cols = np.arange(50+i*nstep, 50+(i+1)*nstep)
            y = []
            x = []
            for col in cols:
                yi = A.image[llim[col]:ulim[col], col]
                y.append(yi / yi.sum())
                x.append(np.arange(llim[col], ulim[col]) - fiber.trace[col])
            x = np.hstack(np.array(x))
            xs = np.sort(x)
            y = np.hstack(np.array(y))
            ax[fi, i].scatter(x, y, alpha=0.1, s=5)
            ax[fi, i].plot(xs, power_law(xs, 0.0004), 'r-')
            ax[fi, i].set_ylim([0.0001, 0.4])
            ax[fi, i].set_yscale('log')
            ax[fi, i].set_xlim([-fsize, fsize])
            if i == 0:
                ax[fi, i].set_ylabel('Fiber %i around y=%i'
                                     % (fib+1, int(np.mean(fiber.trace))))
    f.text(0.5, 0.04, 'Pixels from Fiber Trace', ha='center')
    f.text(0.04, 0.5, 'Normalized Amplitude', va='center', rotation='vertical')
    plt.savefig(outname)
    plt.close(f)
    return mastermaskflat


def relative_throughput(args):
    pass


def write_to_TEX(f, args, overscan, gain, readnoise, darkcounts):
    A = []
    for amp in AMPS:
        A.append(amp)
        A.append(overscan[amp])
        A.append(gain[amp])
        A.append(gain[amp]*readnoise[amp])
    B = []
    for amp in AMPS:
        B.append(amp)
        B.append(darkcounts[amp])
        B.append(darkcounts[amp]*gain[amp])
        B.append(darkcounts[amp]*gain[amp]*600.)
    CreateTex.writeObsSummary(f, A, B)

    obs = ['Bias', 'Darks', 'Pixel flats', 'Photon Transfer Curve']
    mastername = ['masterbias', 'masterdark', 'pixelflat', 'ptc']
    for i, v in enumerate(obs):
        CreateTex.writeImageSummary(f, v)
        A = []
        A.append('%s.png' % (mastername[i]))
        A.append(v)
        CreateTex.writeFigure(f, A)
    obs = ['Masked Flats', 'Fiber Profiles', 'LDLS Spectra']
    mastername = ['mask', 'contrast', 'ldls_spectra']
    for i, v in enumerate(obs):
        CreateTex.writeImageSummary(f, v)
        for amp in AMPS:
            name = mastername[i] + ('_%s.png' % amp)
            A = [name, v + (': %s' % amp)]
            CreateTex.writeFigure(f, A)


def rectify(wave, spec, dl=1., flagv=np.nan):
    wv = np.arange(wave.min(), wave.max(), dl)
    specr = np.zeros((spec.shape[0], len(wv)))
    for i in np.arange(spec.shape[0]):
        specr[i, :] = np.interp(wv, wave[i], spec[i], left=flagv, right=flagv)
    return wv, specr


def main():
    # Read the arguments from the command line
    args = parse_args()
    args = read_in_raw(args)
    # Define output folder
    folder = op.join(args.output, 'CAM_' + args.specid)
    mkpath(folder)

    # Get the bias jumps/structure for each amp
    if not args.dont_check_bias:
        (biasjump_left, biasjump_right, structure,
         overscan, masterbias) = {}, {}, {}, {}, {}
        for amp in AMPS:
            (biasjump_left[amp], biasjump_right[amp],
             structure[amp], overscan[amp],
             masterbias[amp]) = check_bias(args, amp, folder)
        make_plot(masterbias, op.join(folder, 'masterbias.png'))

    # Get the dark jumps/structure and average counts
    if not (args.dont_check_dark or args.dont_check_bias):
        darkcounts, masterdark = {}, {}
        for amp in AMPS:
            darkcounts[amp], masterdark[amp] = check_darks(args, amp, folder,
                                                           masterbias[amp])
        make_plot(masterdark, op.join(folder, 'masterdark.png'))
    else:
        darkcounts, masterdark = {}, {}
        for amp in AMPS:
            darkcounts[amp], masterdark[amp] = (0., 0.)

    # Get the readnoise for each amp
    if not args.dont_check_readnoise:
        readnoise = {}
        for amp in AMPS:
            readnoise[amp] = measure_readnoise(args, amp)

    # Get the gain for each amp
    if not args.dont_check_gain:
        gain = {}
        mn_d = {}
        vr_d = {}
        rd = {}
        for a in AMPS:
            gain[a], mn_d[a], vr_d[a], rd[a] = measure_gain(args, a,
                                                            readnoise[a],
                                                            flow=1,
                                                            fhigh=60000)
        make_ptc_plot(mn_d, vr_d, gain, rd, op.join(folder, 'ptc.png'),
                      lowlim=1, highlim=60000)
    # Get the pixel flat for each amp
    if not args.dont_check_pixelflat:
        masterflat, pixelflat = {}, {}
        for amp in AMPS:
            masterflat[amp], pixelflat[amp] = make_pixelflats(args, amp,
                                                              folder)

        make_plot(pixelflat, op.join(folder, 'pixelflat.png'), vmin=0.95,
                  vmax=1.05)

    if not (args.dont_check_dark or args.dont_check_bias or
            args.dont_check_mask):
        mastermaskflat = {}
        for amp in AMPS:
            mastermaskflat[amp] = check_masked_fibers(args, amp,
                                                      masterbias[amp],
                                                      masterdark[amp],
                                                      op.join(folder,
                                                              'mask_%s.png'
                                                              % amp),
                                                      folder)

    if not (args.dont_check_dark or args.dont_check_bias or
            args.dont_check_ldls):
        masterflat = {}
        for amp in AMPS:
            masterflat[amp] = check_ldls(args, amp, masterbias[amp],
                                         masterdark[amp],
                                         op.join(folder,
                                                 'contrast_%s.png' % amp),
                                         folder, gain[amp])
    # Writing everything to a ".tex" file
    if not (args.dont_check_bias or args.dont_check_dark or
            args.dont_check_readnoise or args.dont_check_gain or
            args.dont_check_pixelflat):
        filename = op.join(folder, 'calibration.tex')
        with open(filename, 'w') as f:
            CreateTex.writeHeader(f, args.specid)
            write_to_TEX(f, args, overscan, gain, readnoise, darkcounts)
            CreateTex.writeEnding(f)


if __name__ == '__main__':
    main()
