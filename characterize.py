#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""
import matplotlib
matplotlib.use('agg')
import argparse as ap
import numpy as np
import textwrap
import glob
import os.path as op
from amplifier import Amplifier
from utils import biweight_location, biweight_midvariance, biweight_filter2d
from CreateTexWriteup import CreateTex
from distutils.dir_util import mkpath
from astropy.io import fits
import matplotlib.pyplot as plt
from operator import itemgetter 
import logging
from scipy.signal import medfilt2d


cmap = plt.get_cmap('Greys_r')

AMPS = ["LL", "LU","RU", "RL"]


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
    description = textwrap.dedent('''Characterize - 
    
                     This script does ... (Fill in Later)
                     
                     ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("--ifuslot", nargs='?', type=str, 
                        help='''Single ifuslot value. [REQUIRED]
                        Ex: "075".''', default = None)

    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''Single specid value. [REQUIRED]
                        Ex: "304".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "camra"
                        Ex: "camra" for lab data,
                            "virus" for mountain.''', default = "camra")

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory
                        Default: \"characterized"''', 
                        default="characterized")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")

    parser.add_argument("-bd","--biadir_date", nargs='?', type=str,
                        help='''Bias Directory Date.    
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-bo","--biadir_obsid", nargs='?', type=str,
                        help='''Bias Directory ObsID.    
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-be","--biadir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-xd","--pxfdir_date", nargs='?', type=str,
                        help='''Pixel Flat Directory Date.    
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-xo","--pxfdir_obsid", nargs='?', type=str,
                        help='''Pixel Flat ObsID.    
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-xe","--pxfdir_expnum", nargs='?', type=str,
                        help='''Pixel Flat exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-pd","--ptcdir_date", nargs='?', type=str,
                        help='''Photon Transfer Curve Directory Date.    
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-po","--ptcdir_obsid", nargs='?', type=str,
                        help=''' Photon Transfer Curve ObsID.    
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-pe","--ptcdir_expnum", nargs='?', type=str,
                        help='''Photon Tranfer Curve exposure number.
                        Ex: \"1,3\" or \"05\"''', default=None) 

    parser.add_argument("-fd","--fltdir_date", nargs='?', type=str,
                        help='''Fiber Flat Directory Date.    
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-fo","--fltdir_obsid", nargs='?', type=str,
                        help='''Fiber Flat ObsID.    
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-fe","--fltdir_expnum", nargs='?', type=str,
                        help='''Fiber Flat exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-dd","--drkdir_date", nargs='?', type=str,
                        help='''Dark Directory Date.   
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-do","--drkdir_obsid", nargs='?', type=str,
                        help='''Dark Directory ObsID.  
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-de","--drkdir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-v","--verbose", help='''Verbose.''',
                        action="count", default=0)

    parser.add_argument("-q","--quick", help='''Quicker Version.''',
                        action="count", default=0)

    parser.add_argument("-dcb","--dont_check_bias", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)

    parser.add_argument("-dcd","--dont_check_dark", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)

    parser.add_argument("-dcr","--dont_check_readnoise", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)
                        
    parser.add_argument("-dcg","--dont_check_gain", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)                        

    parser.add_argument("-dcp","--dont_check_pixelflat", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)

                          
    args = parser.parse_args(args=argv)
    
    return args


def read_in_raw(args):
    log = logging.getLogger('characterize')
    # Check that the arguments are filled
    if args.ifuslot:
        args.ifuslot = "%03d"%int(args.ifuslot)
    else:
        msg = 'No IFUSLOT was provided.'
        log.error(msg)   

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations = []
    if not args.dont_check_bias:
        observations.append('bia')
    if not args.dont_check_dark:
        observations.append('drk')
    if not args.dont_check_gain:
        observations.append('ptc')
    if not args.dont_check_pixelflat:
        observations.append('pxf')
    for obs in observations:
        amp_list = []
        for label in labels[:2]:
            getattr(args, obs+label)
            if getattr(args, obs+label) is None:
                msg = '%s%s was not provided' %(obs, label)
                log.error(msg) 
            else:
                setattr(args, obs+label, 
                        getattr(args, obs+label).replace(" ", "").split(','))
        if getattr(args, obs+labels[2]) is not None:
            setattr(args, obs+labels[2], 
                    getattr(args, obs+labels[2]).replace(" ", "").split(','))
        for date in getattr(args, obs+labels[0]):
            for obsid in getattr(args, obs+labels[1]):
                if getattr(args, obs+labels[2]) is not None:   
                    for expnum in getattr(args, obs+labels[2]):
                        folder = op.join(date, 
                                         args.instr, 
                                         "{:s}{:07d}".format(args.instr,
                                         int(obsid)), 
                                         "exp{:02d}".format(int(expnum)),
                                         args.instr)
                        files = sorted(glob.glob(op.join(args.rootdir, folder, 
                                                         '*_%s*' %args.ifuslot)))
                        for fn in files:
                            amp_list.append(Amplifier(fn, '', name=obs))
                            amp_list[-1].subtract_overscan()
                            amp_list[-1].trim_image()
                            args.specid = amp_list[-1].specid
                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr, 
                                                         int(obsid)))
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*', 
                                                     args.instr, '*_%s*' %args.ifuslot)))
                    for fn in files:
                        amp_list.append(Amplifier(fn, '', name=obs))
                        amp_list[-1].subtract_overscan()
                        amp_list[-1].trim_image()
                        args.specid = amp_list[-1].specid
        setattr(args,obs+'_list', amp_list)

    return args       


def make_plot(image_dict, outfile_name, vmin=-5, vmax=15):
    a,b = image_dict[AMPS[0]].shape
    fig = plt.figure(figsize=((1.*b/a)*4,4))
    for i,amp in enumerate(AMPS):
        ax = plt.subplot(2, 2, i+1)
        ax.imshow(image_dict[amp], vmin=vmin, vmax=vmax, cmap=cmap, 
                  origin='lower', interpolation='none')
        ax.text(b*.1, a*.7, amp, fontsize=24, color='r')
        ax.set_xticks([])
        ax.set_yticks([])
    plt.subplots_adjust(wspace=0.025, hspace=0.025)    
    fig.savefig(outfile_name)

def check_bias(args, amp, folder, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    left_edge, right_edge, structure, overscan = [], [], [], []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    log = args.bia_list[sel[0]].log
    
    overscan_list = [[v.overscan_value for i,v in enumerate(args.bia_list) 
                                   if v.amp == amp]]
    overscan = biweight_location(overscan_list)
    log.info('Overscan value for %s: %0.3f' %(amp, overscan))
    # Loop through the bias list and measure the jump/structure
    big_array = np.array([v.image for v in itemgetter(*sel)(args.bia_list)])
    if args.quick:
        func = np.median
    else:
        func = biweight_location
    masterbias = func(big_array, axis=(0,))
    if not args.quick:
        masterbias = biweight_filter2d(masterbias, (25,5),(5,1))
    #else:
    #    masterbias = medfilt2d(masterbias, [25,5])
    a,b = masterbias.shape
    #masterbias = biweight_filter2d(masterbias, (25,5), (3,1))
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))
    
    log.info('Writing masterbias_%s.fits' %(amp))
    hdu.writeto(op.join(folder, 'masterbias_%s_%s.fits' %(args.specid, amp)), clobber=True)

    left_edge = func(masterbias[:,edge:edge+width])
    right_edge = func(masterbias[:,(b-width-edge):(b-edge)])
    structure = func(masterbias[:,edge:(b-edge)],axis=(0,))
    
    log.info('Left edge - Overscan, Right edge - Overscan: %0.3f, %0.3f'
             %(left_edge, right_edge))
    return left_edge, right_edge, structure, overscan, masterbias


def check_darks(args, amp, folder, masterbias, edge=3, width=10):
    
    # Create empty lists for the left edge jump, right edge jump, and structure
    dark_counts = []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.drk_list) if v.amp == amp]
    log = args.drk_list[sel[0]].log

    if len(sel)<=2 or args.quick:
        func = np.median
    else:
        func = biweight_location
    log.info('Writing masterdark_%s.fits' %(amp))
    big_array = np.array([v.image - masterbias 
                                    for v in itemgetter(*sel)(args.drk_list)])
    masterdark = func(big_array,  axis=(0,))
    a,b = masterdark.shape
    hdu = fits.PrimaryHDU(np.array(masterdark, dtype='float32'))
    hdu.writeto(op.join(folder, 'masterdark_%s_%s.fits' %(args.specid, amp)), clobber=True)

    # Loop through the bias list and measure the jump/structure
    for am in itemgetter(*sel)(args.drk_list):
        a,b = am.image.shape
        dark_counts.append(func(am.image - masterbias) / am.exptime)
    s = biweight_location(dark_counts)
    log.info('Average Dark counts/s: %0.5f' %s)
    return s, masterdark  
    
    
def measure_readnoise(args, amp):
    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    log = args.bia_list[sel[0]].log
    # Make array of all bias images for given amp
    array_images = np.array([bia.image for bia in itemgetter(*sel)(args.bia_list)])
    
    # Measure the biweight midvariance (sigma) for a given pixel and take
    # the biweight average over all sigma to reduce the noise in the first 
    # measurement.
    if args.quick:
        func1 = np.median
        func2 = np.std
    else:
        func1 = biweight_location
        func2 = biweight_midvariance
    S = func1(func2(array_images,axis=(0,)))
    log.info("RDNOISE(ADU) for %s: %01.3f" %(amp, S)) 
    
    return S
    

def measure_gain(args, amp, rdnoise, flow=500, fhigh=35000, fnum=35):
    sel = [i for i,v in enumerate(args.ptc_list) if v.amp == amp]
    log = args.ptc_list[sel[0]].log
    s_sel = list(np.array(sel)[
                 np.array([args.ptc_list[i].basename for i in sel]).argsort()])
    npairs = len(sel) / 2
    a,b = args.ptc_list[sel[0]].image.shape
    array_avg = np.zeros((npairs, a, b))
    array_diff = np.zeros((npairs, a, b))
    if args.quick:
        func1 = np.median
        func2 = np.std
    else:
        func1 = biweight_location
        func2 = biweight_midvariance
    for i in xrange(npairs):
        F1 = args.ptc_list[s_sel[2*i]].image
        F2 = args.ptc_list[s_sel[2*i+1]].image
        m1 = func1(F1)
        m2 = func1(F2)
        array_avg[i,:,:] = (F1 + F2) / 2.
        array_diff[i,:,:] = F1*m2/m1 - F2
    bins = np.logspace(np.log10(flow), np.log10(fhigh), fnum)
    gn = []
    array_avg = array_avg.ravel()
    array_diff = array_diff.ravel()

    for i in xrange(len(bins)-1):
        loc = np.where((array_avg>bins[i]) * (array_avg<bins[i+1]))[0]
        std = func2(array_diff[loc])
        vr   = (std**2 - 2.*rdnoise**2) / 2.
        mn = func1(array_avg[loc])
        log.info("%s | Gain: %01.3f | RDNOISE (e-): %01.3f | <ADU>: %0.1f | "
                  "VAR: %0.1f | Pixels: %i" 
                  %(amp, mn / vr, mn / vr * rdnoise, mn, vr, len(loc))) 
        gn.append(mn/vr)
    s = func1(gn)
    log.info("Average Gain measurement for %s: %0.3f" 
             %(amp, s))
    return s


def make_pixelflats(args, amp, folder):
    sel = [i for i,v in enumerate(args.pxf_list) if v.amp == amp]
    log = args.pxf_list[sel[0]].log

    a,b = args.pxf_list[sel[0]].image.shape
    masterflat = np.zeros((len(sel), a, b))
    
    if args.quick or True:
        for i,am in enumerate(itemgetter(*sel)(args.pxf_list)):
            masterflat[i,:,:] = am.image
        masterflat = np.median(masterflat, axis=(0,))
        smooth = medfilt2d(masterflat, (151,1))
        masterflat = np.where(masterflat==0, 0.0, smooth / masterflat)
        smooth = medfilt2d(masterflat, (1,151))
        pixflat = np.where(masterflat==0, 0.0, smooth / masterflat)
    else:
        print('nothing here')
        #pixflat = np.zeros((len(sel), a, b))
        #for i,am in enumerate(itemgetter(*sel)(args.pxf_list)):

        #masterflat = biweight_filter2d(masterflat, (151,11), (1,1))
        
        
        #log.info('Dividing masterflat for %s' %amp)
        #    pixflat[i,:,:] = np.where(masterflat==0, 0.0, am.image / masterflat)
        #pixflat = func(pixflat, axis=(0,))
        #pixflat = 
        #a,b = pixflat.shape
    hdu = fits.PrimaryHDU(np.array(pixflat, dtype='float32'))
    log.info('Writing pixelflat_%s.fits' %amp)
    hdu.writeto(op.join(folder, 'pixelflat_%s.fits' %amp), clobber=True)

    return masterflat, pixflat

def relative_throughput(args):
    pass

def write_to_TEX(f, args, overscan, gain, readnoise, darkcounts):
    A = []
    for amp in AMPS:
        A.append(amp)
        A.append(overscan[amp])
        A.append(gain[amp])
        A.append(gain[amp]*readnoise[amp])
    B = []
    for amp in AMPS:
        B.append(amp)
        B.append(darkcounts[amp])
        B.append(darkcounts[amp]*gain[amp])
        B.append(darkcounts[amp]*gain[amp]*600.)
    CreateTex.writeObsSummary(f, A, B)
    
    obs = ['Bias', 'Darks', 'Pixel flats']
    mastername = ['masterbias', 'masterdark', 'pixelflat']
    for i, v in enumerate(obs):
        A = [v]
        A.append('%s.png' %(mastername[i]))
        A.append(v)
        CreateTex.writeImageSummary(f, A)
    
def main():
    # Read the arguments from the command line
    args = parse_args()    
    args = read_in_raw(args)
    # Define output folder
    folder = op.join(args.output,'CAM_'+args.specid)
    mkpath(folder)
    
    # Get the bias jumps/structure for each amp
    if not args.dont_check_bias:
        (biasjump_left, biasjump_right, structure, 
         overscan, masterbias) = {}, {}, {}, {}, {}
        for amp in AMPS:
            (biasjump_left[amp], biasjump_right[amp], 
             structure[amp], overscan[amp], 
             masterbias[amp]) = check_bias(args, amp, folder)
        make_plot(masterbias, op.join(folder, 'masterbias.png'))
    
    # Get the dark jumps/structure and average counts
    if not (args.dont_check_dark or args.dont_check_bias):
        darkcounts, masterdark = {}, {}
        for amp in AMPS:
            darkcounts[amp], masterdark[amp] = check_darks(args, amp, folder,
                                                           masterbias[amp])
        make_plot(masterdark, op.join(folder, 'masterdark.png'))
    
    # Get the readnoise for each amp
    if not args.dont_check_readnoise:
        readnoise = {}
        for amp in AMPS:
            readnoise[amp] = measure_readnoise(args, amp)
    
    # Get the gain for each amp
    if not args.dont_check_gain:
        gain = {}
        for amp in AMPS:
            gain[amp] = measure_gain(args, amp, readnoise[amp])

    # Get the pixel flat for each amp
    if not args.dont_check_pixelflat:
        masterflat, pixelflat= {}, {}
        for amp in AMPS:
            masterflat[amp], pixelflat[amp] = make_pixelflats(args, amp, folder)
            
        make_plot(pixelflat, op.join(folder, 'pixelflat.png'), vmin=0.95, vmax=1.05)
    # Writing everything to a ".tex" file
    if not (args.dont_check_bias or args.dont_check_dark 
            or args.dont_check_readnoise or args.dont_check_gain
            or args.dont_check_pixelflat):
        filename = op.join(folder,'calibration.tex')
        with open(filename,'w') as f:
            CreateTex.writeHeader(f, args.specid)
            write_to_TEX(f, args, overscan, gain, readnoise, darkcounts)
            CreateTex.writeEnding(f)
        
if __name__ == '__main__':
    main()  