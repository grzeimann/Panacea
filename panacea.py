#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IFU Reduction Code 
------------------
Built for the VIRUS instrument as well as LRS2 on HET

Fibers

"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import textwrap
import glob
from distutils.dir_util import mkpath
import time

import argparse as ap
import numpy as np
import pandas as pd
import os.path as op
from astropy.io import fits

from amplifier import Amplifier
from fiber_utils import get_model_image

Amps = ["LL", "RU"]
Amp_dict = {"LL": ["LU","L"], "RU": ["RL","R"]}

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
    description = textwrap.dedent('''Panacea - 
    
                     This script does ....
                     
                     ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "020,008".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data.''', default = "virus")

    parser.add_argument("-d","--debug", help='''Debug.''',
                        action="count", default=0)

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory
                        Default: \"virus_reductions"''', 
                        default="virus_reductions")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")

    parser.add_argument("--configdir", nargs='?', type=str, 
                        help='''Config Directory
                        Default: \"/work/03946/hetdex/maverick/virus_config\"''', 
                        default="/work/03946/hetdex/maverick/virus_config")

    parser.add_argument("--biasdir", type=str,
                        help= "Directory of biases to use",
                        default="/work/03946/hetdex/maverick/virus_config/lib_bias")
 
    parser.add_argument("--darkdir", type=str,
                        help= "Directory of darks to use",
                        default="/work/03946/hetdex/maverick/virus_config/lib_dark")
                       
    parser.add_argument("-sd","--scidir_date", nargs='?', type=str,
                        help='''Science Directory Date.     [REQUIRED]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                        help='''Science Directory ObsID.    [REQUIRED]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-td","--twidir_date", nargs='?', type=str,
                        help='''Twi Directory Date.     [REQUIRED]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-to","--twidir_obsid", nargs='?', type=str,
                        help='''Twi Directory ObsID.    [REQUIRED]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-te","--twidir_expnum", nargs='?', type=str,
                        help='''Twi Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        msg = 'No SPECID was provided.'
        parser.error(msg)   

    if args.scidir_date is None:
        msg = 'No science directory date was provided'
        parser.error(msg) 
    else:
        args.scidir_date = args.scidir_date.replace(" ", "").split(',')

    if args.scidir_obsid is None:
        msg = 'No science directory ObsID was provided'
        parser.error(msg) 
    else:
        args.scidir_obsid = args.scidir_obsid.replace(" ", "").split(',')

    if args.scidir_expnum is not None:
        args.scidir_expnum = args.zrodir_expnum.replace(" ", "").split(',')
    
    args.sci_df = pd.DataFrame(columns=['Files', 'Output', 'Amp', 'Specid', 
                                        'Ifuslot', 'Ifuid'])     
    cnt=0
    for date in args.scidir_date:
        for obsid in args.scidir_obsid:
            if args.scidir_expnum is not None:   
                for expnum in args.scidir_expnum:
                    folder = op.join(date, 
                                    args.instr, 
                                    "{:s}{:07d}".format(args.instr,int(obsid)), 
                                    "exp{:02d}".format(int(expnum)),
                                    args.instr)
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*')))
                    if files:
                        mkpath(op.join(args.output,folder))   
                    for fn in files:
                        F = fits.open(fn)
                        outfolder = op.join(args.output,folder)    
                        amp = (F[0].header['CCDPOS'].replace(' ', '') 
                                         + F[0].header['CCDHALF'].replace(' ', ''))
                        sp = '%03d' %F[0].header['SPECID']
                        ifuid = F[0].header['IFUID'].replace(' ', '')
                        ifuslot = '%03d' %F[0].header['IFUSLOT']
                        args.sci_df.loc[cnt] = pd.Series({'Files':fn, 
                                                          'Output':outfolder, 
                                                          'Specid':sp,
                                                          'Ifuslot': ifuslot,
                                                          'Ifuid': ifuid, 
                                                          'Amp': amp})
                        cnt+=1
            else:
                folder = op.join(date, args.instr,
                                 "{:s}{:07d}".format(args.instr,int(obsid)))
                files = sorted(glob.glob(op.join(args.rootdir, folder, '*', args.instr, 
                                                 '*')))
                if files:
                    nfiles = sorted(glob.glob(op.join(args.output, folder, '*')))
                    for nfile in nfiles:
                        mkpath(op.join(nfile, args.instr))
                for fn in files:
                    F = fits.open(fn)
                    exp = op.basename(op.dirname(op.dirname(fn)))
                    outfolder = op.join(args.output,folder, exp, args.instr)
                    amp = (F[0].header['CCDPOS'].replace(' ', '') 
                                     + F[0].header['CCDHALF'].replace(' ', ''))
                    sp = '%03d' %F[0].header['SPECID']
                    ifuid = F[0].header['IFUID'].replace(' ', '')
                    ifuslot = '%03d' %F[0].header['IFUSLOT']
                    args.sci_df.loc[cnt] = pd.Series({'Files':fn, 
                                                      'Output':outfolder, 
                                                      'Specid':sp,
                                                      'Ifuslot': ifuslot,
                                                      'Ifuid': ifuid, 
                                                      'Amp': amp})
                    cnt+=1 

    if args.twidir_date is None:
        msg = 'No twi directory date was provided'
        parser.error(msg) 
    else:
        args.twidir_date = args.twidir_date.replace(" ", "").split(',')

    if args.twidir_obsid is None:
        msg = 'No twi directory ObsID was provided'
        parser.error(msg) 
    else:
        args.twidir_obsid = args.twidir_obsid.replace(" ", "").split(',')

    if args.twidir_expnum is not None:
        args.twidir_expnum = args.zrodir_expnum.replace(" ", "").split(',')
        
        
    args.twi_df = pd.DataFrame(columns=['Files', 'Output', 'Amp', 'Specid', 
                                        'Ifuslot', 'Ifuid']) 
    cnt=0
    for date in args.twidir_date:
        for obsid in args.twidir_obsid:
            if args.twidir_expnum is not None:   
                for expnum in args.twidir_expnum:
                    folder = op.join(date, 
                                    args.instr, 
                                    "{:s}{:07d}".format(args.instr,int(obsid)), 
                                    "exp{:02d}".format(int(expnum)),
                                    args.instr)
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*')))
                    if files:
                        mkpath(op.join(args.output,folder))                        
                    for fn in files:
                        F = fits.open(fn)
                        outfolder = op.join(args.output,folder)
                        amp = (F[0].header['CCDPOS'].replace(' ', '') 
                                         + F[0].header['CCDHALF'].replace(' ', ''))
                        sp = '%03d' %F[0].header['SPECID']
                        ifuid = F[0].header['IFUID'].replace(' ', '')
                        ifuslot = '%03d' %F[0].header['IFUSLOT']
                        args.twi_df.loc[cnt] = pd.Series({'Files':fn,
                                                          'Output':outfolder, 
                                                          'Specid':sp,
                                                          'Ifuslot': ifuslot,
                                                          'Ifuid': ifuid, 
                                                          'Amp': amp})
                        cnt+=1 
            else:
                folder = op.join(date, args.instr,
                                 "{:s}{:07d}".format(args.instr,int(obsid)))
                files = sorted(glob.glob(op.join(args.rootdir, folder, '*', args.instr, 
                                                 '*')))
                if files:
                    nfiles = sorted(glob.glob(op.join(args.output, folder, '*')))
                    for nfile in nfiles:
                        mkpath(op.join(nfile, args.instr))
                for fn in files:
                    F = fits.open(fn)
                    exp = op.basename(op.dirname(op.dirname(fn)))
                    outfolder = op.join(args.output,folder, exp, args.instr)
                    amp = (F[0].header['CCDPOS'].replace(' ', '') 
                                     + F[0].header['CCDHALF'].replace(' ', ''))
                    sp = '%03d' %F[0].header['SPECID']
                    ifuid = F[0].header['IFUID'].replace(' ', '')
                    ifuslot = '%03d' %F[0].header['IFUSLOT']
                    args.twi_df.loc[cnt] = pd.Series({'Files':fn, 
                                                      'Output':outfolder,
                                                      'Specid':sp,
                                                      'Ifuslot': ifuslot,
                                                      'Ifuid': ifuid, 
                                                      'Amp': amp})
                    cnt+=1
    return args 

def main():
    args = parse_args()
    if args.debug:
        t1 = time.time()
    for spec in args.specid:
        spec_ind_twi = np.where(args.twi_df['Specid'] == spec)[0]
        spec_ind_sci = np.where(args.sci_df['Specid'] == spec)[0]
        # Calibration
        for amp in Amps:
            amp_ind_twi = np.where(args.twi_df['Amp'] == amp)[0]
            twi_sel = np.intersect1d(spec_ind_twi, amp_ind_twi)
            amp_ind_sci = np.where(args.sci_df['Amp'] == amp)[0]
            sci_sel = np.intersect1d(spec_ind_sci, amp_ind_sci)
            for ind in twi_sel:
                twi1 = Amplifier(args.twi_df['Files'][ind],
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, darkpath=args.biasdir,
                                 virusconfig=args.configdir)
                twi1.load_fibers()
                if twi1.fibers:
                    if twi1.fibers[0].fiber_to_fiber is not None:
                        print("Loading Cal for %s, %s: %s" %(spec, amp, 
                                                    args.twi_df['Files'][ind]))
                else:
                    twi1.get_fiber_to_fiber(use_default_profile=False, 
                               init_lims=[3490,5500], interactive=False,
                               check_fibermodel=True, check_wave=True)
                twi2 = Amplifier(args.twi_df['Files'][ind].replace(amp, Amp_dict[amp][0]),
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, darkpath=args.biasdir,
                                 virusconfig=args.configdir)
                twi2.load_fibers()
                if twi2.fibers:
                    if twi2.fibers[0].fiber_to_fiber is not None:
                        print("Loading Cal for %s, %s: %s" %(spec, Amp_dict[amp][0], 
                                                    args.twi_df['Files'][ind]))
                else:
                    twi2.get_fiber_to_fiber(use_default_profile=False, 
                               init_lims=[3490,5500], interactive=False,
                               check_fibermodel=True, check_wave=True)
                image1 = get_model_image(twi1.image, twi1.fibers, 'fiber_to_fiber',
                                        debug=twi1.debug)
                image2 = get_model_image(twi1.image, twi1.fibers, 'fiber_to_fiber',
                                        debug=twi1.debug)
                a,b = image1.shape
                new = np.zeros((a*2,b))
                new[:a,:] = image1
                new[a:,:] = image2
                hdu = fits.PrimaryHDU(new, header=twi1.header)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'mastertrace_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    Amp_dict[amp][1]))
                hdu.writeto(outname)
                twi1.save_fibers()
                twi2.save_fibers()
    if args.debug:
        t2=time.time()
        print("Total Time taken: %0.2f s" %(t2-t1))

    
        
    


if __name__ == '__main__':
    main()    



        

       