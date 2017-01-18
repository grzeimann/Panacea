# -*- coding: utf-8 -*-
"""
Argument Parser
------------------
Built for the VIRUS instrument as well as LRS2 on HET

@author: gregz
"""

import pandas as pd
import argparse as ap
from distutils.dir_util import mkpath
import textwrap
import glob
import sys
import os.path as op
from astropy.io import fits

import config

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
    
                     This script does ... (Fill in Later)
                     
                     ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("-rt","--reduce_twi", 
                        help='''Reduce Twighlight frames for calibration''',
                        action="count", default=0)

    parser.add_argument("-rs","--reduce_sci", 
                        help='''Reduce Science frames''',
                        action="count", default=0)
                        
    parser.add_argument("-uos","--use_other_sky", 
                        help='''Use another frame for the sky_spectrum''',
                        action="count", default=0)

    parser.add_argument("-at","--adjust_trace", 
                        help='''Adjust trace using science frame.''',
                        action="count", default=0)
    
    parser.add_argument("-rff","--refit_fiber_to_fiber", 
                        help='''Using the sky spectrum, refit the fiber to
                                fiber.''',
                        action="count", default=0)

    parser.add_argument("-p","--pixelflats", 
                        help='''Turn off using pixel flats.''',
                        action="count", default=0)

    parser.add_argument("-sfs","--start_from_scratch", 
                        help='''Re-fiberextract, sky-subtract, cosmic ray reject.''',
                        action="count", default=0)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. [REQUIRED]
                        Ex: "020,008".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data,
                            "lrs2" for lrs2.''', default = "virus")
                            
    parser.add_argument("--instr_side", nargs='?', type=str, 
                        help='''Instrument side to process. 
                        Default: "blue"
                        Ex: "blue" for LRS2B,
                            "red" for LRS2R.''', default = "blue")

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory
                        Default: \"reductions"''', 
                        default="reductions")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")

    parser.add_argument("--configdir", nargs='?', type=str, 
                        help='''Config Directory
                        Default: \"/work/03946/hetdex/maverick/virus_config\"''', 
                        default="/work/03946/hetdex/maverick/virus_config")

    parser.add_argument("--biasdir", type=str,
                        help='''Bias Library
                        Default: \"/work/03946/hetdex/maverick/virus_config/lib_bias\"''', 
                        default="/work/03946/hetdex/maverick/virus_config/lib_bias")
 
    parser.add_argument("--darkdir", type=str,
                        help='''Dark Library
                        Default: \"/work/03946/hetdex/maverick/virus_config/lib_dark\"''',
                        default="/work/03946/hetdex/maverick/virus_config/lib_dark")
                       
    parser.add_argument("-sd","--scidir_date", nargs='?', type=str,
                        help='''Science Directory Date.     [REQUIRED, if --reduce_sci]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                        help='''Science Directory ObsID.    [REQUIRED, if --reduce_sci]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-td","--twidir_date", nargs='?', type=str,
                        help='''Twi Directory Date.     [REQUIRED, if --reduce_twi]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-to","--twidir_obsid", nargs='?', type=str,
                        help='''Twi Directory ObsID.    [REQUIRED, if --reduce_twi]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-te","--twidir_expnum", nargs='?', type=str,
                        help='''Twi Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-skd","--skydir_date", nargs='?', type=str,
                        help='''Sky Directory Date.     [REQUIRED, --use_other_sky]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-sko","--skydir_obsid", nargs='?', type=str,
                        help='''Sky Directory ObsID.    [REQUIRED, --use_other_sky]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-ske","--skydir_expnum", nargs='?', type=str,
                        help='''Sky Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-d","--debug", help='''Debug.''',
                        action="count", default=0)
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        msg = 'No SPECID was provided.'
        parser.error(msg)   
    
    if args.instr.lower() == 'virus':
        args.wvl_dict = config.virus_wl
        args.specname = config.virus_sn
        args.fsize = config.virus_fs
        args.disp = config.virus_di
        args.ifucen_fn = config.virus_fn
        args.cube_scale = config.virus_cs
        args.fibmodel_bins = config.virus_bn
        args.fibmodel_sig = config.virus_sig
        args.fibmodel_pow = config.virus_pow
        args.dark_mult = config.virus_dm
    if args.instr.lower() == 'lrs2':
        if args.instr_side.lower() == 'blue':
            args.wvl_dict = config.lrs2b_wl
            args.specname = config.lrs2b_sn
            args.fsize = config.lrs2b_fs
            args.disp = config.lrs2b_di
            args.ifucen_fn = config.lrs2b_fn
            args.cube_scale = config.lrs2b_cs
            args.fibmodel_bins = config.lrs2b_bn
            args.fibmodel_sig = config.lrs2b_sig
            args.fibmodel_pow = config.lrs2b_pow
            args.dark_mult = config.lrs2b_dm
        else:
            args.wvl_dict = config.lrs2r_wl
            args.specname = config.lrs2r_sn
            args.fsize = config.lrs2r_fs
            args.disp = config.lrs2r_di    
            args.ifucen_fn = config.lrs2r_fn
            args.cube_scale = config.lrs2r_cs
            args.fibmodel_bins = config.lrs2r_bn
            args.fibmodel_sig = config.lrs2r_sig
            args.fibmodel_pow = config.lrs2r_pow
            args.dark_mult = config.lrs2r_dm



    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations=[]
    args.check_if_twi_exists=False
    if args.reduce_sci:
        observations.append('sci')
        if not args.reduce_twi:
            observations.append('twi')
            args.check_if_twi_exists = True
    if args.reduce_twi:
        observations.append('twi')
    if args.use_other_sky:
        observations.append('sky')
    for obs in observations:
        for label in labels[:2]:
            getattr(args, obs+label)
            if getattr(args, obs+label) is None:
                msg = '%s%s was not provided' %(obs, label)
                parser.error(msg) 
            else:
                setattr(args, obs+label, 
                        getattr(args, obs+label).replace(" ", "").split(','))
        if getattr(args, obs+labels[2]) is not None:
            setattr(args, obs+labels[2], 
                    getattr(args, obs+labels[2]).replace(" ", "").split(','))
        DF =  pd.DataFrame(columns=['Files', 'Output', 'Amp', 'Specid', 
                                    'Ifuslot', 'Ifuid'])     
        cnt=0
        for date in getattr(args, obs+labels[0]):
            for obsid in getattr(args, obs+labels[1]):
                if getattr(args, obs+labels[2]) is not None:   
                    for expnum in getattr(args, obs+labels[2]):
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
                            DF.loc[cnt] = pd.Series({'Files':fn, 
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
                        DF.loc[cnt] = pd.Series({'Files':fn, 'Output':outfolder, 
                                                 'Specid':sp, 'Ifuslot': ifuslot,
                                                 'Ifuid': ifuid, 'Amp': amp})
                        cnt+=1
        setattr(args, obs+'_df', DF)
        
    if args.reduce_sci:
        cals=['twi']
        if args.use_other_sky:
            cals.append('sky')
        else: 
            args.sky_dir = None
        for cal in cals:
            if getattr(args, cal+labels[0]) is None:
                print("Please provide one "+cal+labels[0])
                sys.exit(1) 
            if len(getattr(args, cal+labels[0]))>1:
                print("Please provide only one "+cal+labels[0])
                print("I am cowardly quitting instead of making a smart program.")
                sys.exit(1)
            if getattr(args, cal+labels[1]) is None:
                print("Please provide one "+cal+labels[1])
                sys.exit(1) 
            if len(getattr(args, cal+labels[1]))>1:
                print("Please provide only one "+cal+labels[1])
                print("I am cowardly quitting instead of making a smart program.")
                sys.exit(1)
            if getattr(args, cal+labels[2]) is None:
                print("Please provide one "+cal+labels[2])
                sys.exit(1) 
            if len(getattr(args, cal+labels[2]))>1:
                print("Please provide only one "+cal+labels[2])
                print("I am cowardly quitting instead of making a smart program.")
                sys.exit(1)             
            for date in getattr(args, cal+labels[0]):
                for obsid in getattr(args, cal+labels[1]):
                    if getattr(args, cal+labels[2]) is not None:   
                        for expnum in getattr(args, cal+labels[2]):
                            setattr(args, cal+'_dir', op.join(args.output, date, args.instr, 
                                                   "{:s}{:07d}".format(args.instr,int(obsid)), 
                                                   "exp{:02d}".format(int(expnum)),
                                                   args.instr))
            
        

    return args