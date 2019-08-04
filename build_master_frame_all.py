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
from math_utils import biweight

karl_tarlist = '/work/00115/gebhardt/maverick/gettar/%starlist'

def write_fits(hdu, name):
    try:
        hdu.writeto(name, overwrite=True)
    except:
        hdu.writeto(name, clobber=True)


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
        for daten in args.daterange:
            daten = '%04d%02d%02d' % (daten.year, daten.month, daten.day)
            process = subprocess.Popen('cat %s | grep %s | grep _%sLL | '
                                       'grep %s' %
                                       (tname, args.kind, ifuslot, daten),
                                       stdout=subprocess.PIPE, shell=True)
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                b = line.rstrip()
                c = op.join(args.rootdir, b)
                filenames.append(c[:-14])
    return filenames

def get_tarfiles(filenames):
    tarnames = []
    for fn in filenames:
        tarbase = op.dirname(op.dirname(op.dirname(fn))) + '.tar'
        if tarbase not in tarnames:
            tarnames.append(tarbase)
    return tarnames

def get_ifuslots(tarfolder):
    T = tarfile.open(tarfolder, 'r')
    flag = True
    ifuslots = []
    expnames = []
    while flag:
        a = T.next()
        try:
            name = a.name
        except:
            flag = False
        if name[-5:] == '.fits':
            namelist = name.split('/')
            expn = namelist[-3]
            ifuslot = namelist[-1].split('_')[1][:3]
            if ifuslot not in ifuslots:
                ifuslots.append(ifuslot)
            if expn not in expnames:
                expnames.append(expn)
            if len(expnames) > 1:
                flag = False
    return sorted(ifuslots)

def get_unique_ifuslots(tarfolders):
    dates = [tarfolder.split('/')[-2] for tarfolder in tarfolders]
    utars, inds = np.unique(dates, return_index=True)
    utarfolders = [tarfolders[ind] for ind in inds]
    ifuslots = []
    for tarfolder in utarfolders:
        result = get_ifuslots(tarfolder)
        ifuslots = list(np.union1d(ifuslots, result))
    return ifuslots

def get_image(fn):
    tarbase = op.dirname(op.dirname(op.dirname(fn))) + '.tar'
    if op.exists(tarbase):
        T = tarfile.open(tarbase, 'r')
        s = '/'.join(fn.split('/')[-4:])
        fn = T.extractfile(s)
    
    A = Amplifier(fn, '')
    ly = A.biassec[2]
    hy = A.biassec[3]
    lx = A.biassec[0]+1
    hx = A.biassec[1]
    A.overscan_value = biweight_location(A.image[ly:hy, lx:hx])
    A.image[:] = A.image - A.overscan_value
    A.trim_image()
    tlist = A.header['UT'].split(':')
    hour, minute, second = (int(tlist[0]), int(tlist[1]),
                            int(tlist[2].split('.')[0]))
    d = datetime(A.date.year, A.date.month, A.date.day, hour, minute, second)
    return (A.image * 1., A.specid, '%04d%02d%02d' % (A.date.year,
                                                      A.date.month, 
                                                      A.date.day),
            A.header['OBJECT'], A.header, d)


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
    if args.kind == 'cmp':
        mname = 'mastercmp'
        objnames = ['hg', 'cd-a']
    sname = 'masterstd'
    bia_list = []
    for itm in file_list:
        fn = itm + '%s%s_%s.fits' % (ifuslot, amp, args.kind)
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
    if args.kind != 'cmp':
        big_array = np.array([v[0] for v in bia_list])
        if args.kind == 'flt':
                norm = np.median(big_array, axis=(1, 2))
                sel = norm > 0.
                big_array = big_array[sel] / norm[sel, np.newaxis, np.newaxis]
                big_array *= np.median(norm)
        func = biweight
        masterbias, masterstd = func(big_array, nan_treatment=False, calc_std=True,
                                     axis=(0,))
    else:
        masterbias = np.zeros(bia_list[0][0].shape)
        masterstd = np.zeros(bia_list[0][0].shape)

        for objname in objnames:
            big_array = np.array([v[0] for v in bia_list
                                  if v[3].lower() == objname])
            func = biweight
            masterim, masterst = func(big_array, nan_treatment=False, calc_std=True,
                                         axis=(0,))
            masterbias += masterim
            masterstd += masterst / len(objnames)
    for masterim, Name in zip([masterbias, masterstd], [mname, sname]):
        a, b = masterim.shape
        hdu = fits.PrimaryHDU(np.array(masterim, dtype='float32'),
                              header=bia_list[0][4])
        
        d1 = bia_list[0][5]
        d2 = bia_list[-1][5]
        d4 = (d1 + timedelta(seconds=(d2-d1).seconds / 2.) +
              timedelta(days=(d2-d1).days / 2.))
        avgdate = '%04d%02d%02dT%02d%02d%02d' % (d4.year, d4.month, d4.day,
                                                 d4.hour, d4.minute, d4.second)
        mkpath(op.join(args.folder, avgdate))
        args.log.info('Writing %s_%s_%s.fits' % (Name, bia_list[-1][1], amp))
        hdu.header['OBJECT'] = '%s-%s' % (bia_list[0][2], bia_list[-1][2])
        write_fits(hdu, op.join(args.folder, avgdate, '%s_%s_%s.fits' %
                   (Name, bia_list[-1][1], amp)))

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
if args.kind not in ['zro', 'drk', 'sci', 'twi', 'flt', 'cmp']:
    args.log.error('"--kind" argument did not match "zro" or "drk"')
    sys.exit(1)

filenames = get_filenames(args)
tarnames = get_tarfiles(filenames)
ifuslots = get_unique_ifuslots(tarnames)

args.log.info('Number of unique ifuslots: %i' % len(ifuslots))

for ifuslot in ifuslots:
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
            args.log.warning('No files found for %s on %s' % (ifuslot, date))
            break
        chunks = np.array_split(filenames, nbins)
        for chunk in chunks:
            datestr = op.basename(chunk[0])[:8]
            build_master_frame(chunk, ifuslot, amp, args, datestr)
