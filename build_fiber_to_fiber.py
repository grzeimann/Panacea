# -*- coding: utf-8 -*-
"""
Build new Fiber to Fiber
"""

import glob
import numpy as np
import os.path as op

from astropy.io import fits
from distutils.dir_util import mkpath
from fiber_utils import bspline_x0
from input_utils import setup_parser, set_daterange, setup_logging
from scipy.interpolate import splev, splrep
from utils import biweight_location


def check_if_type(date, obsid, args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation and exposure.
    '''
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
    Build directory structure and search for unique observations, and return
    a single file for each observation and exposure.
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
            F = fits.open(name)
            for i, attribute in enumerate(attributes):
                s[i].append(F[attribute].data)
        except IOError:
            args.log.warning('%s not found, filling with zeros' % name)
            for i, attribute in enumerate(attributes):
                s[i].append(np.zeros((112, 1032)))
    return [np.vstack(si) for si in s]


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

    ifu_spline_dict = {}
    # HARDCODED SIZE FOR SPEED BUT MUST MATCH SIZE OF "rw" BELOW.
    B, c = bspline_x0(np.linspace(0, 1, 2001), nknots=7)
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
            for ifu in ifus:
                if ifu not in ifu_spline_dict:
                    ifu_spline_dict[ifu] = []
            for exposure in exps:
                file_list = [fn for fn, e in zip(filenames, e_list)
                             if e == exposure]
                ifuslot_list = [i for i, e in zip(i_list, e_list)
                                if e == exposure]
                args.log.info('Building Fiber to Fiber for %s, observation %s,'
                              ' exposure %s' % (date, obsid, exposure))
                allspec = ([])
                for filen, ifu in zip(file_list, ifuslot_list):
                    wave, spec = grab_attribute(filen, args,
                                                attributes=['wavelength',
                                                            'spectrum'])
                    rw, rs = rectify(wave, spec, minwave=3500., maxwave=5500.)
                    allspec.append(rs)
                rectwave = rw
                allspec = np.array(allspec)
                avgspec = np.nanmedian(allspec, axis=(0, 1))
                for sp, ifu in zip(allspec, ifuslot_list):
                    div = sp / avgspec
                    splinecoeff = np.zeros((sp.shape[0], c.shape[1]))
                    for i, d in enumerate(div):
                        splinecoeff[i, :] = np.linalg.lstsq(c, d,
                                                            rcond=None)[0]
                    ifu_spline_dict[ifu].append(splinecoeff)

    ifu_ftf_dict = {}
    mkpath(args.outdir)
    for ifu in ifu_spline_dict:
        fibers_spline_coeff = biweight_location(np.array(ifu_spline_dict[ifu]),
                                                axis=(0,))
        ftf = np.zeros(rs.shape)
        for i, fiber in enumerate(fibers_spline_coeff):
            ftf[i, :] = np.dot(c, fiber)
        ifu_ftf_dict[ifu] = ftf
        P = fits.PrimaryHDU(np.array(ftf, dtype='float32'))
        P.writeto(args.outdir, '%s_ftf.fits' % ifu, overwrite=True)

if __name__ == '__main__':
    main()
