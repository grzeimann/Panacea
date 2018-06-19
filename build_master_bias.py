# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:32:27 2018

@author: gregz
"""

import glob
import os.path as op
import numpy as np

from amplifier import Amplifier
from astropy.io import fits
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from utils import biweight_location


def write_fits(hdu, name):
    try:
        hdu.writeto(name, overwrite=True)
    except:
        hdu.writeto(name, clobber=True)


def build_filenames(date, args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    basedir = op.join(args.rootdir, date, args.instrument,
                      args.instrument + '0000*', 'exp*', args.instrument)
    filenames = sorted(glob.glob(op.join(basedir, '2*_zro.fits')))
    dirnames = [op.dirname(fn) for fn in filenames]
    unique_dirnames, ind = np.unique(dirnames, return_index=True)
    return [filenames[i][:-14] for i in ind]


def build_master_frame(file_list, ifuslot, amp, args, date):
    # Create empty lists for the left edge jump, right edge jump, and structure

    bia_list = []
    for itm in file_list:
        fn = itm + '%s%s_zro.fits' % (ifuslot, amp)
        bia_list.append(Amplifier(fn, ''))
        bia_list[-1].subtract_overscan()
        bia_list[-1].trim_image()

    # Select only the bias frames that match the input amp, e.g., "RU"
    if not len(bia_list):
        args.log.warning('No bias frames found for date range given')
        return None

    # Loop through the bias list and measure the jump/structure
    big_array = np.array([v.image for v in bia_list])
    func = biweight_location
    masterbias = func(big_array, axis=(0,))

    a, b = masterbias.shape
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))
    mkpath(op.join(args.folder, date))
    args.log.info('Writing masterbias_%s_%s.fits' % (bia_list[-1].specid, amp))
    write_fits(hdu, op.join(args.folder, date, 'masterbias_%s_%s.fits' %
               (bia_list[-1].specid, amp)))

parser = setup_parser()
parser.add_argument("-f", "--folder",
                    help='''Output folder''',
                    type=str, default='masterbias')

args = parser.parse_args(args=None)
args.log = setup_logging(logname='build_master_bias')
args = set_daterange(args)
filenames = []
for date in args.daterange:
    date = '%04d%02d%02d' % (date.year, date.month, date.day)
    filenames = filenames + build_filenames(date, args)
for ifuslot in ['056', '066']:
    for amp in ['LL', 'LU', 'RL', 'RU']:
        build_master_frame(filenames, ifuslot, amp, args,
                           '%04d%02d%02d' % (args.daterange[0].year,
                                             args.daterange[0].month,
                                             args.daterange[0].day))
