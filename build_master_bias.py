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


def get_image(fn):
    A = Amplifier(fn, '')
    A.subtract_overscan()
    A.trim_image()
    return A.image * 1., A.specid


def build_master_frame(file_list, ifuslot, amp, args, date, maxnum):
    # Create empty lists for the left edge jump, right edge jump, and structure

    bia_list = []
    for itm in file_list:
        fn = itm + '%s%s_zro.fits' % (ifuslot, amp)
        try:
            bia_list.append(get_image(fn))
        except:
            args.log.warning('Could not load %s' % fn)

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

parser.add_argument("-m", "--maxnum",
                    help='''Maximum number of bias frames in masterbias''',
                    type=int, default=400)

args = parser.parse_args(args=None)
args.log = setup_logging(logname='build_master_bias')
args = set_daterange(args)
filenames = []
for date in args.daterange:
    date = '%04d%02d%02d' % (date.year, date.month, date.day)
    filenames = filenames + build_filenames(date, args)
for ifuslot in ['056', '066']:
    for amp in ['LL', 'LU', 'RL', 'RU']:
        chunks = np.array_split(filenames,
                                len(filenames) / (args.maxnum+1) + 1)
        for chunk in chunks:
            datestr = op.basename(chunk[0])[:8]
            build_master_frame(chunk, ifuslot, amp, args, datestr)
