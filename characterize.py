#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""
import matplotlib
import argparse as ap
import numpy as np
import glob
import os.path as op
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

matplotlib.use('agg')

matplotlib.rcParams['font.sans-serif'] = "Meiryo"
matplotlib.rcParams['font.family'] = "sans-serif"
# plt.style.use('seaborn-colorblind')

cmap = plt.get_cmap('Greys_r')

AMPS = ["LL", "LU", "RU", "RL"]


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
    obstype = ['bia', 'drk', 'pxf', 'ptc', 'flt']
    obsletter = ['b', 'd', 'x', 'p', 'f']
    obsname = ['Bias', 'Dark', 'Pixel Flat', 'Photon Transfer Curve',
               'Masked Fiber Flat']
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
        args.fltdir_date = args.date
        args.fltdir_obsid = '%03d%04d' % (int(args.specid), 3)

    if not args.dont_check_bias:
        observations.append('bia')
    if not args.dont_check_dark:
        observations.append('drk')
    if not args.dont_check_gain:
        observations.append('ptc')
    if not args.dont_check_pixelflat:
        observations.append('pxf')
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
                        files = sorted(glob.glob(op.join(args.rootdir, folder,
                                                         '*_%s*'
                                                         % args.ifuslot)))
                        for fn in files:
                            amp_list.append(Amplifier(fn, '', name=obs))
                            amp_list[-1].subtract_overscan()
                            amp_list[-1].trim_image()
                            args.specid = amp_list[-1].specid
                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr,
                                                         int(obsid)))
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*',
                                                     args.instr,
                                                     '*_%s*' % args.ifuslot)))
                    for fn in files:
                        amp_list.append(Amplifier(fn, '', name=obs))
                        amp_list[-1].subtract_overscan()
                        amp_list[-1].trim_image()
                        args.specid = amp_list[-1].specid
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

    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i, v in enumerate(args.bia_list) if v.amp == amp]
    log = args.bia_list[sel[0]].log

    overscan_list = [[v.overscan_value for i, v in enumerate(args.bia_list)
                      if v.amp == amp]]
    overscan = biweight_location(overscan_list)
    log.info('Overscan value for %s: %0.3f' % (amp, overscan))
    # Loop through the bias list and measure the jump/structure
    big_array = np.array([v.image for v in itemgetter(*sel)(args.bia_list)])
    if args.quick:
        func = np.median
    else:
        func = biweight_location
    masterbias = func(big_array, axis=(0,))

    a, b = masterbias.shape
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))

    log.info('Writing masterbias_%s.fits' % (amp))
    hdu.writeto(op.join(folder, 'masterbias_%s_%s.fits' % (args.specid, amp)),
                clobber=True)

    left_edge = func(masterbias[:, edge:edge+width])
    right_edge = func(masterbias[:, (b-width-edge):(b-edge)])
    structure = func(masterbias[:, edge:(b-edge)], axis=(0,))

    log.info('Left edge - Overscan, Right edge - Overscan: %0.3f, %0.3f'
             % (left_edge, right_edge))
    return left_edge, right_edge, structure, overscan, masterbias


def check_darks(args, amp, folder, masterbias, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    dark_counts = []

    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i, v in enumerate(args.drk_list) if v.amp == amp]
    log = args.drk_list[sel[0]].log

    if len(sel) <= 2 or args.quick:
        func = np.median
    else:
        func = biweight_location
    log.info('Writing masterdark_%s.fits' % (amp))
    big_array = np.array([v.image - masterbias
                          for v in itemgetter(*sel)(args.drk_list)])
    masterdark = func(big_array,  axis=(0,))
    a, b = masterdark.shape
    hdu = fits.PrimaryHDU(np.array(masterdark, dtype='float32'))
    hdu.writeto(op.join(folder, 'masterdark_%s_%s.fits' % (args.specid, amp)),
                clobber=True)

    # Loop through the bias list and measure the jump/structure
    for am in itemgetter(*sel)(args.drk_list):
        a, b = am.image.shape
        dark_counts.append(func(am.image - masterbias) / am.exptime)
    s = biweight_location(dark_counts)
    log.info('Average Dark counts/s: %0.5f' % s)
    return s, masterdark


def measure_readnoise(args, amp):
    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i, v in enumerate(args.bia_list) if v.amp == amp]
    log = args.bia_list[sel[0]].log
    # Make array of all bias images for given amp
    array_images = np.array([bia.image for bia in
                             itemgetter(*sel)(args.bia_list)])

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
    sel = [i for i, v in enumerate(args.ptc_list) if v.amp == amp]
    log = args.ptc_list[sel[0]].log
    s_sel = list(np.array(sel)[
                 np.array([args.ptc_list[i].basename for i in sel]).argsort()])
    npairs = len(sel) / 2
    a, b = args.ptc_list[sel[0]].image.shape
    array_avg = np.zeros((npairs, a, b))
    array_diff = np.zeros((npairs, a, b))
    if args.quick:
        func1 = np.median
        func2 = np.std
    else:
        func1 = biweight_location
        func2 = biweight_midvariance
    for i in xrange(npairs):
        F1 = args.ptc_list[s_sel[2*i]].image
        F2 = args.ptc_list[s_sel[2*i+1]].image
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
        std = func2(array_diff[loc])
        vr = (std**2) / 2.
        vr_c = (std**2 - 2.*rdnoise**2) / 2.
        mn = func1(array_avg[loc])
        log.info("%s | Gain: %01.3f | RDNOISE (e-): %01.3f | <ADU>: %0.1f | "
                 "VAR: %0.1f | Pixels: %i"
                 % (amp, mn / vr_c, mn / vr_c * rdnoise, mn, vr, len(loc)))
        gn.append(mn / vr_c)
        mn_list.append(mn)
        vr_list.append(vr)
    sel = np.where((np.array(mn_list) > 1000.)*(np.array(mn_list) < 15000.))[0]
    s = func1(np.array(gn)[sel])
    log.info("Average Gain measurement for %s: %0.3f"
             % (amp, s))
    return s, mn_list, vr_list, rdnoise**2


def make_pixelflats(args, amp, folder):
    sel = [i for i, v in enumerate(args.pxf_list) if v.amp == amp]
    log = args.pxf_list[sel[0]].log

    a, b = args.pxf_list[sel[0]].image.shape
    masterflat = np.zeros((len(sel), a, b))

    for i, am in enumerate(itemgetter(*sel)(args.pxf_list)):
        masterflat[i, :, :] = am.image
    masterflat = np.median(masterflat, axis=(0,))
    smooth = medfilt2d(masterflat, (151, 1))
    masterflat = np.where(masterflat < 1e-8, 0.0, smooth / masterflat)
    smooth = medfilt2d(masterflat, (1, 151))
    pixflat = np.where(masterflat < 1e-8, 0.0, smooth / masterflat)

    hdu = fits.PrimaryHDU(np.array(pixflat, dtype='float32'))
    log.info('Writing pixelflat_%s.fits' % amp)
    hdu.writeto(op.join(folder, 'pixelflat_%s.fits' % amp), clobber=True)

    return masterflat, pixflat


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
        A = [v]
        A.append('%s.png' % (mastername[i]))
        A.append(v)
        CreateTex.writeImageSummary(f, A)


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
