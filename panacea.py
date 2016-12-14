#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IFU Reduction Code 
------------------
Built for the VIRUS instrument as well as LRS2 on HET

Incomplete Documentation



"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
                        
import matplotlib
matplotlib.use('agg')
import textwrap
import glob
from distutils.dir_util import mkpath
import time

import sys
import argparse as ap
import numpy as np
import pandas as pd
import os.path as op
from astropy.io import fits
from pyhetdex.cure.distortion import Distortion


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
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. [REQUIRED]
                        Ex: "020,008".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data,
                            "lrs2" for lrs2.''', default = "virus")

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

    parser.add_argument("-d","--debug", help='''Debug.''',
                        action="count", default=0)
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        msg = 'No SPECID was provided.'
        parser.error(msg)   

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations=[]
    if args.reduce_sci:
        observations.append('sci')
        if not args.reduce_twi:
            observations.append('twi')
    if args.reduce_twi:
        observations.append('twi')
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
        if getattr(args, 'twi'+labels[0]) is None:
            print("Please provide one "+"twi"+labels[0])
            sys.exit(1) 
        if len(getattr(args, 'twi'+labels[0]))>1:
            print("Please provide only one "+"twi"+labels[0])
            print("I am cowardly quitting instead of making a smart program.")
            sys.exit(1)
        if getattr(args, 'twi'+labels[1]) is None:
            print("Please provide one "+"twi"+labels[1])
            sys.exit(1) 
        if len(getattr(args, 'twi'+labels[1]))>1:
            print("Please provide only one "+"twi"+labels[1])
            print("I am cowardly quitting instead of making a smart program.")
            sys.exit(1)
        if getattr(args, 'twi'+labels[2]) is None:
            print("Please provide one "+"twi"+labels[2])
            sys.exit(1) 
        if len(getattr(args, 'twi'+labels[2]))>1:
            print("Please provide only one "+"twi"+labels[2])
            print("I am cowardly quitting instead of making a smart program.")
            sys.exit(1)             
        for date in getattr(args, 'twi'+labels[0]):
            for obsid in getattr(args, 'twi'+labels[1]):
                if getattr(args, 'twi'+labels[2]) is not None:   
                    for expnum in getattr(args, obs+labels[2]):
                        args.cal_dir = op.join(args.output, date, args.instr, 
                                               "{:s}{:07d}".format(args.instr,int(obsid)), 
                                               "exp{:02d}".format(int(expnum)),
                                               args.instr)
                        
                else:
                    print("You need to provide an exposure number for twi if "
                          "you are reducing science frames.")
                    print("I am cowardly quitting instead of making a smart "
                          "program.")
                    sys.exit(1)
    return args

def matrixCheby2D_7(x, y):
    if isinstance(x, (tuple, list)):
        x = np.asarray(x)
    if isinstance(y, (tuple, list)):
        y = np.asarray(y)

    T2x = 2. * x**2 - 1.
    T3x = 4. * x**3 - 3. * x
    T4x = 8. * x**4 - 8. * x**2 + 1.
    T5x = 16. * x**5 - 20. * x**3 + 5. * x
    T6x = 32. * x**6 - 48. * x**4 + 18. * x**2 - 1.
    T7x = 64. * x**7 - 112. * x**5 + 56. * x**3 - 7. * x
    T2y = 2. * y**2 - 1.
    T3y = 4. * y**3 - 3. * y
    T4y = 8. * y**4 - 8. * y**2 + 1.
    T5y = 16. * y**5 - 20. * y**3 + 5. * y
    T6y = 32. * y**6 - 48. * y**4 + 18. * y**2 - 1
    T7y = 64. * y**7 - 112. * y**5 + 56. * y**3 - 7 * y
    
    return np.vstack((T7x, T6x, T5x, T4x, T3x, T2x, x, T7y, T6y, T5y, 
                      T4y, T3y, T2y, y, T6x*y, x*T6y, T5x*T2y, T2x*T5y,
                      T4x*T3y, T3x*T4y, T5x*y, x*T5y, T4x*T2y, T2x*T4y, 
                      T3x*T3y, T4x*y, x*T4y, T3x*T2y, T2x*T3y, T3x*y, 
                      x*T3y, T2x*T2y, T2x*y, x*T2y, x*y, np.ones(x.shape))).swapaxes(0,1)
                      
                      
def recreate_fiberextract(instr1, instr2, wavelim, disp):
    col = int(instr1.D / 2)
    intv = [1, 1+instr1.D]
    ypos = np.array([fiber.trace+intv[v] for v,instr in enumerate([instr1,instr2]) 
                                   for fiber in instr.fibers])
    allspec = np.array([fiber.spectrum/fiber.fiber_to_fiber 
                                   for v,instr in enumerate([instr1,instr2]) 
                                   for fiber in instr.fibers])                                       
    allskys = np.array([(fiber.spectrum-fiber.sky_spectrum)/fiber.fiber_to_fiber 
                                   for v,instr in enumerate([instr1,instr2]) 
                                   for fiber in instr.fibers])
    allwave = np.array([fiber.wavelength for v,instr in enumerate([instr1,instr2]) 
                                   for fiber in instr.fibers])
    f1 = ypos[:,col]
    order = np.argsort(f1)[::-1]
    orderspec = allspec[order,:]
    orderskys = allskys[order,:]
    orderwave = allwave[order,:]
    wv = np.arange(wavelim[0], wavelim[1]+disp, disp)
    a,b = orderspec.shape
    newspec = np.zeros((a, len(wv)))
    newskys = np.zeros((a, len(wv)))
    for i in xrange(a):
        diff_wave = np.diff(orderwave[i,:])
        wi = orderwave[i,:-1]
        df = np.interp(wv, wi, diff_wave, left=0.0, right=0.0)
        fl = np.interp(wv, wi, orderspec[i,:-1], left=0.0, right=0.0)
        fs = np.interp(wv, wi, orderskys[i,:-1], left=0.0, right=0.0)
        newspec[i,:] = np.where(df!=0, fl/df*disp, 0.0)
        newskys[i,:] = np.where(df!=0, fs/df*disp, 0.0)
    return newspec, newskys
        
    
def recalculate_dist_coeff(D, instr1, instr2):
    col = int(instr1.D / 2)
    intv = [1, 1+instr1.D]
    ypos = np.array([fiber.trace+intv[v] for v,instr in enumerate([instr1,instr2]) 
                                   for fiber in instr.fibers])
    xpos = np.array([np.arange(fiber.D) for instr in [instr1,instr2] 
                                        for fiber in instr.fibers])
    f1 = ypos[:,col]
    f0 = np.sort(f1)[::-1]
    fpos = np.zeros((ypos.shape))
    k=0
    for instr in [instr1,instr2]:
        for fiber in instr.fibers:
            fpos[k,:] = f1[k]
            k+=1
            
    wpos = np.array([fiber.wavelength for instr in [instr1,instr2] 
                                      for fiber in instr.fibers])    
    
    Vxy = matrixCheby2D_7(D._scal_x(xpos.ravel()), 
                             D._scal_y(ypos.ravel()))
    Vwf = matrixCheby2D_7(D._scal_w(wpos.ravel()), 
                             D._scal_f(fpos.ravel()))
    Vxf = matrixCheby2D_7(D._scal_x(xpos.ravel()), 
                             D._scal_f(fpos.ravel())) 
    D.reference_f_.data = f0*1.
    D.fiber_par_.data = np.linalg.lstsq(Vxy, fpos.ravel())[0]
    D.wave_par_.data = np.linalg.lstsq(Vxy, wpos.ravel())[0]
    D.x_par_.data = np.linalg.lstsq(Vwf, xpos.ravel())[0]
    D.y_par_.data = np.linalg.lstsq(Vwf, ypos.ravel())[0]
    D.fy_par_.data = np.linalg.lstsq(Vxf, ypos.ravel())[0]
    D.x_offsets = f0*0.
    D.wave_offsets = f0*0.
    return D
    
def make_fiber_image(Fe, header, outname):
    a,b = Fe.shape
    hdu = fits.PrimaryHDU(Fe, header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,a)
    hdu.header['CRVAL1'] = 3500
    hdu.header['CDELT1'] = 1.9
    hdu.header['CD1_1'] = 1.9
    hdu.header['CRPIX1'] = 1
    hdu.writeto(outname, clobber=True)    
    
def make_spectrograph_image(image1, image2, header, outname):
    a,b = image1.shape
    new = np.zeros((a*2,b))
    new[:a,:] = image1
    new[a:,:] = image2
    hdu = fits.PrimaryHDU(new, header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,2*a)
    hdu.writeto(outname, clobber=True)    
            
def reduce_science(args):
    for spec in args.specid:
        spec_ind_sci = np.where(args.sci_df['Specid'] == spec)[0]
        for amp in Amps:
            amp_ind_sci = np.where(args.sci_df['Amp'] == amp)[0]
            sci_sel = np.intersect1d(spec_ind_sci, amp_ind_sci) 
            for ind in sci_sel:
                if args.debug:
                    print("Working on Sci for %s, %s" %(spec, amp))   
                sci1 = Amplifier(args.sci_df['Files'][ind],
                                 args.sci_df['Output'][ind],
                                 calpath=args.cal_dir, 
                                 debug=False, refit=False, dark_mult=1.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir)
                sci1.load_fibers()
                sci1.load_all_cal()
                sci1.sky_subtraction()
                sci2 = Amplifier(args.sci_df['Files'][ind].replace(amp, Amp_dict[amp][0]),
                                 args.sci_df['Output'][ind],
                                 calpath=args.sci_df['Output'][ind], 
                                 debug=False, refit=False, dark_mult=1.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir)
                sci2.load_fibers()
                sci2.load_all_cal()
                print(len(sci2.fibers),args.sci_df['Files'][ind].replace(amp, Amp_dict[amp][0]),
                      args.sci_df['Output'][ind])
                sci2.sky_subtraction()
                outname = op.join(args.sci_df['Output'][ind],
                                  'S%s_%s_sci_%s.fits' %(
                                  op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                  args.sci_df['Ifuslot'][ind], Amp_dict[amp][1]))
                make_spectrograph_image(sci1.clean_image, sci2.clean_image, sci1.header, 
                           outname)
               
                Fe, FeS = recreate_fiberextract(sci1, sci2, wavelim=[3500,5500], 
                                      disp=1.9)
                outname = op.join(args.sci_df['Output'][ind],
                                  'Fe%s_%s_sci_%s.fits' %(
                                  op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                  args.sci_df['Ifuslot'][ind], Amp_dict[amp][1]))
                make_fiber_image(Fe, sci1.header, outname)
                outname = op.join(args.sci_df['Output'][ind],
                                  'FeS%s_%s_sci_%s.fits' %(
                                  op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                  args.sci_df['Ifuslot'][ind], Amp_dict[amp][1]))
                make_fiber_image(FeS, sci1.header, outname)
                sci1.save_fibers()
                sci2.save_fibers()  
                if args.debug:
                    print("Finished working on Sci for %s, %s" %(spec, amp))                    
                

def reduce_twighlight(args):
    D = Distortion(op.join(args.configdir, 'DeformerDefaults', 
                                        'mastertrace_twi_027_L.dist'))   
    for spec in args.specid:
        spec_ind_twi = np.where(args.twi_df['Specid'] == spec)[0]
        for amp in Amps:
            amp_ind_twi = np.where(args.twi_df['Amp'] == amp)[0]
            twi_sel = np.intersect1d(spec_ind_twi, amp_ind_twi)
            for ind in twi_sel:
                if args.debug:
                    print("Working on Cal for %s, %s" %(spec, amp))                    
                twi1 = Amplifier(args.twi_df['Files'][ind],
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=False, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir)
                twi1.load_fibers()
                twi1.get_fiber_to_fiber(use_default_profile=False, 
                               init_lims=[3490,5500], interactive=False,
                               check_fibermodel=True, check_wave=True)
                twi2 = Amplifier(args.twi_df['Files'][ind].replace(amp, Amp_dict[amp][0]),
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=False, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir)
                twi2.load_fibers()
                twi2.get_fiber_to_fiber(use_default_profile=False, 
                               init_lims=[3490,5500], interactive=False,
                               check_fibermodel=True, check_wave=True)
                image1 = get_model_image(twi1.image, twi1.fibers, 'fiber_to_fiber',
                                        debug=twi1.debug)
                image2 = get_model_image(twi2.image, twi2.fibers, 'fiber_to_fiber',
                                        debug=twi1.debug)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'mastertrace_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    Amp_dict[amp][1]))
                make_spectrograph_image(image1, image2, twi1.header, outname)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'mastertwi_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    Amp_dict[amp][1]))  
                make_spectrograph_image(twi1.image, twi2.image, twi1.header, outname)
                D = recalculate_dist_coeff(D, twi1, twi2)
                outname2 = op.join(args.twi_df['Output'][ind], 
                                  'mastertrace_%s_%s.dist' 
                                  %(args.twi_df['Specid'][ind],
                                    Amp_dict[amp][1]))
                D.writeto(outname2)
                twi1.save_fibers()
                twi2.save_fibers()
                if args.debug:
                    print("Finished working on Cal for %s, %s" %(spec, amp))  
                
def main():
    args = parse_args()
    if args.debug:
        t1 = time.time()
    if args.reduce_twi:
        reduce_twighlight(args)
    if args.reduce_sci:
        reduce_science(args)                                        
    if args.debug:
        t2=time.time()
        print("Total Time taken: %0.2f s" %(t2-t1))

if __name__ == '__main__':
    main()    



        

       