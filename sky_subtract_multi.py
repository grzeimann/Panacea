# -*- coding: utf-8 -*-
"""
Sky Subtract
"""

import matplotlib
matplotlib.use('agg')
import glob
import numpy as np
import os.path as op
import fitsio


from astropy.io import fits
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from scipy.interpolate import splev, splrep
from utils import biweight_location



def check_if_type(date, obsid, args):
    ''' Test if header has IMAGETYP '''
    filenames = glob.glob(op.join(args.rootdir, date, args.instrument,
                                  args.instrument + obsid, 'exp01',
                                  args.instrument, 'multi_*_*_*_LL.fits'))
    try:
        kind = fits.open(filenames[0])[0].header['IMAGETYP']
    except:
        args.log.warn('No IMAGETYP in header for %s and observation %s'
                      % (date, obsid))
        return False
    if kind == args.type:
        return True
    else:
        return False


def build_filenames(date, obsid, args):
    '''
    Build directory structure and search for all the files in a given
    observation and exposure.
    '''
    if args.type == 'twi':
        expstr = '01'
    else:
        expstr = '*'
    filenames = glob.glob(op.join(args.rootdir, date, args.instrument,
                                  args.instrument + obsid, 'exp%s' % expstr,
                                  args.instrument, 'multi_*_*_*_LL.fits'))
    ifuslot_list = [op.basename(fn).split('_')[2] for fn in filenames]
    ifuslots = np.unique(ifuslot_list)
    exposure_list = [op.basename(op.dirname(op.dirname(fn)))[3:]
                     for fn in filenames]
    exposures = np.unique(exposure_list)
    return filenames, ifuslots, exposures, ifuslot_list, exposure_list


def grab_attribute(filename, args, attributes=[],
                   amps=['LL', 'LU', 'RU', 'RL']):
    ''' grab specified attributes from multi* file '''
    basename = filename[:-8]
    s = [[] for a in attributes]
    for amp in amps:
        name = basename + '_%s.fits' % amp
        try:
            for i, attribute in enumerate(attributes):
                s[i].append(fitsio.read(name, attribute))
        except IOError:
            args.log.warning('%s not found, filling with zeros' % name)
            for i, attribute in enumerate(attributes):
                s[i].append(np.zeros((112, 1032)))
        for i, attribute in enumerate(attributes):
            if s[i][-1].shape != (112, 1032):
                s[i][-1] = np.zeros((112, 1032))

    return [np.array(si) for si in s]


def rectify(wave, spec, rectified_dlam=1., minwave=None, maxwave=None):
    ''' Rectify spectra to same "rect_wave" '''
    dlam = np.zeros(wave.shape)
    dlam[:, 1:] = np.diff(wave, axis=1)
    dlam[:, 0] = dlam[:, 1]
    if rectified_dlam is None:
        rectified_dlam = np.nanmedian(dlam)

    rect_wave = np.arange(wave.min(), wave.max() + rectified_dlam,
                          rectified_dlam)
    if minwave is not None and maxwave is not None:
        wnew = np.arange(minwave, maxwave + rectified_dlam,
                         rectified_dlam)
    else:
        wnew = rect_wave * 1.
    rect_spec = np.zeros((spec.shape[0], len(wnew)))
    xs = np.linspace(0, 1, len(rect_wave))
    xn = np.interp(wnew, rect_wave, xs)
    for i in np.arange(spec.shape[0]):
        if np.all(spec[i] == 0):
            rect_spec[i, :] = 0.0
        else:
            y = spec[i] / dlam[i]
            xp = np.interp(wave[i], rect_wave, xs)
            tck = splrep(xp, y)
            rect_spec[i, :] = splev(xn, tck)
    rect_wave = wnew * 1.
    return rect_wave, rect_spec


def main():
    parser = setup_parser()
    parser.add_argument("-t", "--type",
                        help='''Observation Type, twi or sci''',
                        type=str, default='twi')
    parser.add_argument("-o", "--outdir",
                        help='''Out directory for fiber to fiber''',
                        type=str, default='temp')
    args = parser.parse_args(args=None)
    args.log = setup_logging(logname='build_ftf')
    args = set_daterange(args)

    # HARDCODED SIZE FOR SPEED BUT MUST MATCH SIZE OF "rw" BELOW.
    for datet in args.daterange:
        date = '%04d%02d%02d' % (datet.year, datet.month, datet.day)
        obsids = glob.glob(op.join(args.rootdir, date, args.instrument,
                                   args.instrument + '*'))
        obsids = [obsid[-7:] for obsid in obsids]
        for obsid in obsids:
            if not check_if_type(date, obsid, args):
                continue
            filenames, ifus, exps, i_list, e_list = build_filenames(date,
                                                                    obsid,
                                                                    args)
            for exposure in exps:
                file_list = [fn for fn, e in zip(filenames, e_list)
                             if e == exposure]
                ifuslot_list = [i for i, e in zip(i_list, e_list)
                                if e == exposure]
                ifuslot_amp = ['%s%s' % (ifu, amp) for ifu in ifuslot_list
                               for amp in ['LL', 'LU', 'RU', 'RL']]

                args.log.info('Building Fiber to Fiber for %s, observation %s,'
                              ' exposure %s' % (date, obsid, exposure))
                allspec, ftf = ([], [])
                for filen, ifu in zip(file_list, ifuslot_list):
                    args.log.info('Reading in %s' % filen)
                    amps = ['LL', 'LU', 'RU', 'RL']
                    wave, spec, FtF = grab_attribute(filen, args, attributes=[
                                                     'wavelength', 'spectrum',
                                                     'fiber_to_fiber'],
                                                     amps=amps)
                    for wv, sp, amp, Ftf in zip(wave, spec, amps, FtF):
                        rw, rs = rectify(wv, sp, minwave=3500., maxwave=5500.)
                        allspec.append(rs)
                        ftf.append(Ftf)
                allspec = np.array(allspec)
                avgspec = np.nanmedian(allspec, axis=(0, 1))
                ifuslot_offset_dict = {}
                for i in np.arange(0, 2000, 20):
                    cols = np.arange(i, i+20)
                    y = np.median(ftf[:, :, cols], axis=2)
                    y2 = np.median(allspec[:, :, cols] / avgspec[cols], axis=2)
                    for j, ifuamp in enumerate(ifuslot_amp):
                        if ifuamp not in ifuslot_offset_dict:
                            ifuslot_offset_dict[ifuamp] = []
                        offset = y2[j, :] - y[j, :]
                        offsetm = np.median(offset)
                        ifuslot_offset_dict[ifuamp].append(offsetm)

if __name__ == '__main__':
    main()
