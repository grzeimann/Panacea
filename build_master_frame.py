# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:32:27 2018

@author: gregz
"""

import glob
import os.path as op
import numpy as np
import subprocess
import sys
import tarfile

from amplifier import Amplifier
from astropy.io import fits
from datetime import datetime, timedelta
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from utils import biweight_location

karl_tarlist = '/work/00115/gebhardt/maverick/gettar/%starlist'

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
    tarfolder = op.join(args.rootdir, date, args.instrument, 
                        "{:s}0000*.tar".format(args.instrument))
    tarfolders = glob.glob(tarfolder)
    if len(tarfolders):
        filenames = []
        for tarfolder in tarfolders:
            args.log.info('Inspecting %s' % tarfolder)
            T = tarfile.open(tarfolder, 'r')
            flag = True
            while flag:
                a = T.next()
                try:
                    name = a.name
                except:
                    flag = False
                if name[-5:] == '.fits':
                    flag = False
                    if name[-9:] == ('_%s.fits' % args.kind):
                        names = sorted(T.getnames())
                        for namei in names:
                            if namei[-5:] == '.fits':
                                if namei[-9:] == ('_%s.fits' % args.kind):
                                    filenames.append(op.join(op.dirname(
                                                     tarfolder), namei))
    else:
        basedir = op.join(args.rootdir, date, args.instrument,
                          args.instrument + '0000*', 'exp*', args.instrument)
        filenames = sorted(glob.glob(op.join(basedir,
                                             '2*_%s.fits' % args.kind)))
    filenames = sorted(filenames)
    dirnames = [op.dirname(fn) for fn in filenames]
    unique_dirnames, ind = np.unique(dirnames, return_index=True)
    return [filenames[i][:-14] for i in ind]

def get_filenames(args):
    filenames = []
    ifuslot = '%03d' % int(args.ifuslot)
    dates = []
    for date in args.daterange:
        date = '%04d%02d' % (date.year, date.month)
        if date not in dates:
            dates.append(date)
    
    for date in dates:
        tname = karl_tarlist % date
        process = subprocess.Popen('cat %s | grep %s | grep _%s' %
                                   (tname, args.kind, ifuslot),
                                   stdout=subprocess.PIPE, shell=True)
        while True:
            line = process.stdout.readline()
            if not line:
                break
            b = line.rstrip()
            c = op.join(args.rootdir, b)
            filenames.append(c)
    return filenames

def get_image(fn):
    tarbase = op.dirname(op.dirname(op.dirname(fn))) + '.tar'
    if op.exists(tarbase):
        T = tarfile.open(tarbase, 'r')
        s = '/'.join(fn.split('/')[-4:])
        fn = T.extractfile(T.getmember(s))
    A = Amplifier(fn, '')
    ly = A.biassec[2]
    hy = A.biassec[3]
    lx = A.biassec[0]+1
    hx = A.biassec[1]
    A.overscan_value = biweight_location(A.image[ly:hy, lx:hx])
    A.image[:] = A.image - A.overscan_value
    A.trim_image()
    return (A.image * 1., A.specid, '%04d%02d%02d' % (A.date.year,
                                                      A.date.month, 
                                                      A.date.day))


def build_master_frame(file_list, ifuslot, amp, args, date):
    # Create empty lists for the left edge jump, right edge jump, and structure
    if args.kind == 'zro':
        mname = 'masterbias'
    if args.kind == 'drk':
        mname = 'masterdark'
    if args.kind == 'sci':
        mname = 'mastersci'
    if args.kind == 'twi':
        mname = 'mastertwi'
    if args.kind == 'flt':
        mname = 'masterflt'
    bia_list = []
    for itm in file_list:
        fn = itm #+ '%s%s_%s.fits' % (ifuslot, amp, args.kind)
        try:
            bia_list.append(get_image(fn))
        except:
            args.log.warning('Could not load %s' % fn)

    # Select only the bias frames that match the input amp, e.g., "RU"
    if not len(bia_list):
        args.log.warning('No bias frames found for date range given')
        return None

    # Loop through the bias list and measure the jump/structure
    uspec = np.unique([v[1] for v in bia_list])
    if len(uspec) > 1:
        args.log.warning('More than one spectrograph for given ifuslot, '
                         'cowardly exiting.')
        args.log.warning(uspec)
        return None
    big_array = np.array([v[0] for v in bia_list])
    func = biweight_location
    masterbias = func(big_array, axis=(0,))

    a, b = masterbias.shape
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))
    
    d1 = datetime(int(bia_list[0][2][:4]), int(bia_list[0][2][4:6]),
                  int(bia_list[0][2][6:]))
    d2 = datetime(int(bia_list[-1][2][:4]), int(bia_list[-1][2][4:6]),
                  int(bia_list[-1][2][6:]))
    d3 = d1 + timedelta(days=(d2-d1).days/2)
    avgdate = '%04d%02d%02d' % (d3.year, d3.month, d3.day)
    mkpath(op.join(args.folder, avgdate))
    args.log.info('Writing %s_%s_%s.fits' % (mname, bia_list[-1][1], amp))
    hdu.header['OBJECT'] = '%s-%s' % (bia_list[0][2], bia_list[-1][2])
    write_fits(hdu, op.join(args.folder, avgdate, '%s_%s_%s.fits' %
               (mname, bia_list[-1][1], amp)))

parser = setup_parser()
parser.add_argument("-f", "--folder", help='''Output folder''', type=str,
                    default='masterbias')

parser.add_argument("-m", "--maxnum",
                    help='''Maximum number of bias frames in masterbias''',
                    type=int, default=100)

parser.add_argument("-k", "--kind", help='''drk or zro''',  type=str,
                    default='zro')

parser.add_argument("-i", "--ifuslot",  help='''IFUSLOT''', type=str,
                    default='056')


args = parser.parse_args(args=None)
args.log = setup_logging(logname='build_master_bias')
args = set_daterange(args)
args.kind = args.kind.lower()
if args.kind not in ['zro', 'drk', 'sci', 'twi', 'flt']:
    args.log.error('"--kind" argument did not match "zro" or "drk"')
    sys.exit(1)
#filenames = []
#for date in args.daterange:
#    date = '%04d%02d%02d' % (date.year, date.month, date.day)
#    filenames = filenames + build_filenames(date, args)

filenames = get_filenames(args)

for amp in ['LL', 'LU', 'RL', 'RU']:
    date = args.daterange[0]
    date = '%04d%02d%02d' % (date.year, date.month, date.day)
    date_end = args.daterange[-1]
    date_end = '%04d%02d%02d' % (date_end.year, date_end.month, date_end.day)
    args.log.info('Length of filenames for %s-%s: %i' %
                  (date, date_end, len(filenames)))
    if (len(filenames) % args.maxnum) == 0:
        nbins = len(filenames) / args.maxnum
    else:
        nbins = len(filenames) / args.maxnum + 1
    if nbins == 0:
        args.log.warning('No files found for %s on %s' % (args.ifuslot, date))
        break
    chunks = np.array_split(filenames, nbins)
    for chunk in chunks:
        datestr = op.basename(chunk[0])[:8]
        build_master_frame(chunk, args.ifuslot, amp, args, datestr)