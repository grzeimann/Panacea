# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:52:35 2017

@author: gregz
"""
import glob
import re

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import fits
from input_utils import setup_parser, set_daterange, setup_logging


def build_path(reduction_folder, instr, date, obsid, expn):
    folder = op.join(date, instr, "{:s}{:07d}".format(instr, int(obsid)),
                     "exp{:02d}".format(int(expn)), instr)
    return op.join(reduction_folder, folder)


def get_files(args):
    args.log.info('Looking at %i days' % len(args.daterange))
    for date in args.daterange:
        datestr = '%04d%02d%02d' % (date.year, date.month, date.day)
        files = glob.glob(op.join(args.rootdir, datestr, args.instrument, '*',
                          'exp*', args.instrument, 'multi_*.fits'))
    return files


class LRS2Fiber(tb.IsDescription):
    obsind = tb.Int32Col()
    fibnum = tb.Int32Col()
    x = tb.Float32Col()
    y = tb.Float32Col()
    spectrum = tb.Float32Col((2064,))
    wavelength = tb.Float32Col((2064,))
    fiber_to_fiber = tb.Float32Col((2064,))
    sky_spectrum = tb.Float32Col((2064,))


class LRS2Amp(tb.IsDescription):
    obsind = tb.Int32Col()
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(4)
    specid = tb.StringCol(3)
    amp = tb.StringCol(2)
    date = tb.StringCol(8)
    mjd = tb.Float32Col()
    obsid = tb.StringCol(7)
    ra = tb.Float32Col()
    dec = tb.Float32Col()
    pa = tb.Float32Col()
    expn = tb.Int32Col()
    time = tb.StringCol(7)
    ambtemp = tb.Float32Col()
    humidity = tb.Float32Col()
    dewpoint = tb.Float32Col()
    pressure = tb.Float32Col()
    exptime = tb.Float32Col()


def append_file_to_table(fib, amp, F, cnt):
    if 'spectrum' in F:
        n = F['spectrum'].data.shape[0]
        d = F['spectrum'].data.shape[1]
    else:
        return False
    attr = ['spectrum', 'wavelength', 'fiber_to_fiber', 'sky_spectrum']
    for i in np.arange(n):
        fib['obsind'] = cnt
        fib['fibnum'] = i
        if 'ifupos' in F:
            fib['x'] = F['ifupos'].data[i, 0]
            fib['y'] = F['ifupos'].data[i, 1]
        else:
            fib['x'] = -999.0
            fib['y'] = -999.0
        for att in attr:
            if att in F:
                fib[att] = F[att].data[i, :]
            else:
                fib[att] = np.zeros((d,))
        fib.append()
    amp['obsind'] = cnt
    amp['ifuslot'] = '%03d' % int(F[0].header['IFUSLOT'])
    amp['ifuid'] = '%03d' % int(F[0].header['IFUID'])
    amp['specid'] = '%03d' % int(F[0].header['SPECID'])
    amp['amp'] = '%s' % F[0].header['amp'][:2]
    amp['date'] = ''.join(F[0].header['DATE-OBS'].split('-'))
    amp['time'] = ''.join(re.split('[:,.]', F[0].header['UT']))[:7]
    amp['mjd'] = F[0].header['MJD']
    amp['obsid'] = '%07d' % F[0].header['OBSID']
    amp['ra'] = F[0].header['TRAJCRA'] * 15.
    amp['dec'] = F[0].header['TRAJCDEC']
    amp['pa'] = F[0].header['PARANGLE']
    amp['ambtemp'] = F[0].header['AMBTEMP']
    amp['humidity'] = F[0].header['HUMIDITY']
    amp['dewpoint'] = F[0].header['DEWPOINT']
    amp['pressure'] = F[0].header['BAROMPRE']
    amp['exptime'] = F[0].header['EXPTIME']
    amp['expn'] = int(op.basename(op.dirname(op.dirname(F.filename())))[-2:])
    amp.append()
    return True


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parent_parser = setup_parser()

    # Use "--append" to append to an existing file "--outfilename"
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               parents=[parent_parser, ])

    parser.add_argument('-o', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)
    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)

    args = parser.parse_args(argv)
    args.log = setup_logging()

    # Get the daterange over which reduced files will be collected
    args = set_daterange(args)
    files = get_files(args)

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, 'a')
        does_exist = True
    else:
        fileh = tb.open_file(args.outfilename, 'w')
        group = fileh.create_group(fileh.root, 'Info',
                                   'LRS2 Fiber Data and Metadata')
        table1 = fileh.create_table(group, 'Fibers', LRS2Fiber, 'Fiber Info')
        table2 = fileh.create_table(group, 'Amps', LRS2Amp, 'Amp Info')

    # Grab the fiber table and amplifier table for writing
    fibtable = fileh.root.Info.Fibers
    amptable = fileh.root.Info.Amps
    if does_exist:
        cnt = amptable[-1]['obsind']
    else:
        cnt = 1

    for fn in files:
        args.log.info('Working on %s' % fn)
        F = fits.open(fn)
        fib = fibtable.row
        amp = amptable.row
        success = append_file_to_table(fib, amp, F, cnt)
        if success:
            cnt += 1
            fibtable.flush()
            amptable.flush()
        attr = ['spectrum', 'wavelength', 'fiber_to_fiber', 'sky_spectrum',
                'ifupos']
        for att in attr:
            if att in F:
                del F[att].data

        F.close()
    fileh.close()


if __name__ == '__main__':
    main()
