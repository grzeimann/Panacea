#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""
import matplotlib
matplotlib.use('TkAgg')
import argparse as ap
import numpy as np
import textwrap
import glob
import os.path as op
from amplifier import Amplifier
from utils import biweight_location, biweight_midvariance, biweight_filter2d
from progressbar import ProgressBar
from CreateTexWriteup import CreateTex
from distutils.dir_util import mkpath
from astropy.io import fits
import matplotlib.pyplot as plt
from operator import itemgetter 

cmap = plt.get_cmap('Greys')

AMPS = ["LL", "LU", "RU", "RL"]


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
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. [REQUIRED]
                        Ex: "020,008".''', default = None)

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

    parser.add_argument("-d","--debug", help='''Debug.''',
                        action="count", default=0)

    parser.add_argument("-dcb","--dont_check_bias", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)

    parser.add_argument("-dcd","--dont_check_dark", 
                        help='''Don't make masterbias.''',
                        action="count", default=0)

                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = "%03d"%int(args.specid)
    else:
        msg = 'No SPECID was provided.'
        parser.error(msg)   

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations=['bia', 'drk', 'pxf', 'ptc', 'flt']
    for obs in observations:
        amp_list = []
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
                                                         '*')))
                        for fn in files:
                            amp_list.append(Amplifier(fn, '', name=obs))
                            amp_list[-1].subtract_overscan()
                            amp_list[-1].trim_image()
                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr, 
                                                         int(obsid)))
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*', 
                                                     args.instr, '*')))
                    for fn in files:
                        amp_list.append(Amplifier(fn, '', name=obs))
                        amp_list[-1].subtract_overscan()
                        amp_list[-1].trim_image()
        setattr(args,obs+'_list', amp_list)

    return args       


def check_bias(args, amp, folder, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    left_edge, right_edge, structure, overscan = [], [], [], []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    overscan = biweight_location([[v.overscan_value 
                                   for i,v in enumerate(args.bia_list) 
                                   if v.amp == amp]])

    # Loop through the bias list and measure the jump/structure
    masterbias = biweight_location(np.array([v.image 
                                    for v in itemgetter(*sel)(args.bia_list)]), 
                                   axis=(0,))
    a,b = masterbias.shape
    #masterbias = biweight_filter2d(masterbias, (25,5), (3,1))
    hdu = fits.PrimaryHDU(np.array(masterbias, dtype='float32'))
    hdu.writeto(op.join(folder, 'masterbias_%s.fits' %amp), clobber=True)
    s = biweight_location(masterbias[:,edge:(b-edge)], axis=(0,))
    vran = np.max(s)-np.min(s)
    vmin = np.min(s)-0.05*vran
    vmax = np.max(s)+0.05*vran
    plt.figure(figsize=(8,6))
    plt.imshow(masterbias, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower',
               interpolation='none')
    plt.text(b*.1, a*.9, amp, fontsize=24, color='r')
    plt.savefig(op.join(folder, 'masterbias_%s.png' %amp), dpi=150)
    plt.close()
    
    for am in itemgetter(*sel)(args.bia_list):
        left_edge.append(biweight_location(am.image[:,edge:edge+width]))
        right_edge.append(biweight_location(
                                         am.image[:,(b-width-edge):(b-edge)]))
        structure.append(biweight_location(am.image[:,edge:(b-edge)],axis=(0,)))
    return left_edge, right_edge, structure, overscan, masterbias


def check_darks(args, amp, folder, masterbias, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    dark_counts = []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.drk_list) if v.amp == amp]
    if len(sel)>2:
        masterdark = biweight_location(np.array([v.image - masterbias 
                                    for v in itemgetter(*sel)(args.drk_list)]), 
                                       axis=(0,))
        a,b = masterdark.shape
        hdu = fits.PrimaryHDU(np.array(masterdark, dtype='float32'))
        hdu.writeto(op.join(folder, 'masterdark_%s.fits' %amp), clobber=True)
        s = biweight_location(masterdark[:,edge:(b-edge)], axis=(0,))
        vran = np.max(s)-np.min(s)
        vmin = np.min(s)-0.05*vran
        vmax = np.max(s)+0.05*vran
        plt.figure(figsize=(8,6))
        plt.imshow(masterdark, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower',
               interpolation='none')
        plt.text(b*.1, a*.9, amp, fontsize=24, color='r')
        plt.savefig(op.join(folder, 'masterdark_%s.png' %amp), dpi=150)
        plt.close()
    else:
        masterdark = None
       
    
    
    # Loop through the bias list and measure the jump/structure
    for am in itemgetter(*sel)(args.drk_list):
        a,b = am.image.shape
        dark_counts.append(biweight_location(am.image - masterbias) / am.exptime)
    return dark_counts, masterdark  
    
    
def measure_readnoise(args, amp):
    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    
    # Make array of all bias images for given amp
    array_images = np.array([bia.image for bia in itemgetter(*sel)(args.bia_list)])
    
    # Measure the biweight midvariance (sigma) for a given pixel and take
    # the biweight average over all sigma to reduce the noise in the first 
    # measurement.
    S = biweight_location(biweight_midvariance(array_images,axis=(0,)))
    
    print("\n%s | RDNOISE (ADU): %01.3f" %(amp, S)) 
    
    return S
    

def measure_gain(args, amp, rdnoise, flow=500, fhigh=35000, fnum=35):
    sel = [i for i,v in enumerate(args.ptc_list) if v.amp == amp]
    s_sel = list(np.array(sel)[
                 np.array([args.ptc_list[i].basename for i in sel]).argsort()])
    npairs = len(sel) / 2
    a,b = args.ptc_list[sel[0]].image.shape
    array_avg = np.zeros((npairs, a, b))
    array_diff = np.zeros((npairs, a, b))
    for i in xrange(npairs):
        F1 = args.ptc_list[s_sel[2*i]].image
        F2 = args.ptc_list[s_sel[2*i+1]].image
        m1 = biweight_location(F1)
        m2 = biweight_location(F2)
        array_avg[i,:,:] = F1 + F2
        array_diff[i,:,:] = F1*m2/m1 - F2
    bins = np.logspace(np.log10(flow), np.log10(fhigh), fnum)
    gn = []
    array_avg = array_avg.ravel()
    array_diff = array_diff.ravel()

    for i in xrange(len(bins)-1):
        loc = np.where((array_avg>bins[i]) * (array_avg<bins[i+1]))[0]
        print(array_diff[loc])
        std = biweight_midvariance(array_diff[loc])
        vr   = (std**2 - 2.*rdnoise**2) / 2.
        mn = biweight_location(array_avg[loc])
        plt.figure(figsize=(8,6))
        plt.hist(array_diff[loc], 50, range=[mn-3*std,mn+3*std])
        plt.show()
        raw_input("Waiting for you to press enter.")
        plt.close()
        print("%s | Gain: %01.3f | RDNOISE (e-): %01.3f | <ADU>: %0.1f | "
              "VAR: %0.1f | Pixels: %i" 
                %(amp, mn / vr, mn / vr * rdnoise, mn, vr, len(loc))) 
        gn.append(mn/vr)
    return biweight_location(gn)


def make_pixelflats(args, amp, folder):
    sel = [i for i,v in enumerate(args.pxf_list) if v.amp == amp]
    a,b = args.pxf_list[sel[0]].image.shape
    masterflat = np.zeros((len(sel), a, b))
    pixflat = np.zeros((len(sel), a, b))
    for am in itemgetter(*sel)(args.pxf_list):
        masterflat[i,:,:] = biweight_filter2d(am.image, (25,5), (3,1))
    masterflat = biweight_location(masterflat, axis=(0,))
    for am in itemgetter(*sel)(args.pxf_list):
        pixflat[i,:,:] = am.image / masterflat
    pixflat = biweight_location(pixflat, axis=(0,))
    
    a,b = pixflat.shape
    hdu = fits.PrimaryHDU(np.array(pixflat, dtype='float32'))
    hdu.writeto(op.join(folder, 'pixelflat_%s.fits' %amp), clobber=True)

    plt.figure(figsize=(8,6))
    plt.imshow(pixflat, vmin=0.95, vmax=1.05, cmap=cmap, origin='lower',
               interpolation='none')
    plt.text(b*.1, a*.9, amp, fontsize=24, color='r')
    plt.savefig(op.join(folder, 'pixelflat_%s.png' %amp), dpi=150)
    plt.close()
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
    
    amps = ["LU","RL","LL","RU"]
    obs = ['Bias', 'Darks', 'Pixel flats']
    mastername = ['masterbias', 'masterdark', 'pixelflat']
    for i, v in enumerate(obs):
        A = [v]
        for amp in amps:
            A.append(amp)
            A.append('%s_%s.png' %(mastername[i], amp))
        CreateTex.writeImageSummary(f, A)
    
def main():
    # Read the arguments from the command line
    args = parse_args()
    
    # Define output folder
    folder = op.join(args.output,'CAM_'+args.specid)
    mkpath(folder)
    
    # Get the bias jumps/structure for each amp
    if not args.dont_check_bias:
        (biasjump_left, biasjump_right, structure, 
         overscan, masterbias) = {}, {}, {}, {}, {}
        progress = ProgressBar(len(AMPS), 'Checking Biases for %s' %args.specid, 
                               fmt=ProgressBar.FULL)
        for amp in AMPS:
            (biasjump_left[amp], biasjump_right[amp], 
             structure[amp], overscan[amp], 
             masterbias[amp]) = check_bias(args, amp, folder)
            progress.current+=1
            progress()
        progress.done()
        
    
    # Get the dark jumps/structure and average counts
    if not args.dont_check_dark:
        darkcounts, masterdark = {}, {}
        progress = ProgressBar(len(AMPS), 'Checking Darks for %s' %args.specid, 
                               fmt=ProgressBar.FULL)
        for amp in AMPS:
            darkcounts[amp], masterdark[amp] = check_darks(args, amp, folder,
                                                           masterbias[amp])
            progress.current+=1
            progress()
        progress.done()
    
    # Get the readnoise for each amp
    readnoise = {}
    progress = ProgressBar(len(AMPS), 'Measuring Readnoise for %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        readnoise[amp] = measure_readnoise(args, amp)
        progress.current+=1
        progress()
    progress.done()
    
    # Get the gain for each amp
    gain = {}
    progress = ProgressBar(len(AMPS), 'Measuring Gain for %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        gain[amp] = measure_gain(args, amp, readnoise[amp])
        progress.current+=1
        progress()
    progress.done()

    # Get the pixel flat for each amp
    masterflat, pixelflat= {}, {}
    progress = ProgressBar(len(AMPS), 'Measuring pixel flat for %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        masterflat[amp], pixelflat[amp] = make_pixelflats(args, amp, 
                                                          readnoise[amp])
        progress.current+=1
        progress()
    progress.done()
    
    # Writing everything to a ".tex" file

    filename = op.join(folder,'calibration.tex')
    with open(filename) as f:
        CreateTex.writeHeader(f)
        write_to_TEX(f, args, overscan, gain, readnoise, darkcounts)
        CreateTex.writeEnding(f)
        
if __name__ == '__main__':
    main()  