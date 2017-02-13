#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""

import argparse as ap
import numpy as np
import textwrap
import glob
import os.path as op
from amplifier import Amplifier
from utils import biweight_location, biweight_midvariance, biweight_filter2d
from progressbar import ProgressBar

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
                        Default: "virus"
                        Ex: "camra" for lab data,
                            "lrs2" for lrs2.''', default = "camra")

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
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
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


def check_bias(args, amp, edge=3, width=10):
    # Create empty lists for the left edge jump, right edge jump, and structure
    left_edge, right_edge, structure, overscan = [], [], [], []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    overscan = biweight_location([[v.overscan_value 
                                   for i,v in enumerate(args.bia_list) 
                                   if v.amp == amp]])

    # Loop through the bias list and measure the jump/structure
    masterbias = biweight_location(np.array([v.image 
                                             for v in args.bia_list[sel]]), 
                                   axis=(0,))
    a,b = masterbias.shape
    masterbias = biweight_filter2d(masterbias[:,i], (5,25), (1,3))
    for am in args.bia_list[sel]:
        left_edge.append(biweight_location(am.image[:,edge:edge+width]))
        right_edge.append(biweight_location(
                                         am.image[:,(b-width-edge):(b-edge)]))
        structure.append(biweight_location(am[:,edge:(b-edge)],axis=(0,)))
    return left_edge, right_edge, structure, overscan, masterbias


def check_darks(args, amp, masterbias):
    # Create empty lists for the left edge jump, right edge jump, and structure
    dark_counts = []
    
    # Select only the bias frames that match the input amp, e.g., "RU"   
    sel = [i for i,v in enumerate(args.drk_list) if v.amp == amp]
    if len(sel)>2:
        masterdark = biweight_location(np.array([v.image - masterbias 
                                                 for v in args.drk_list[sel]]), 
                                       axis=(0,))
    else:
        masterdark = None
    # Loop through the bias list and measure the jump/structure
    for am in args.drk_list[sel]:
        a,b = am.image.shape
        dark_counts.append(biweight_location(am - masterbias))
    return dark_counts, masterdark  
    
    
def measure_readnoise(args, amp):
    # Select only the bias frames that match the input amp, e.g., "RU"
    sel = [i for i,v in enumerate(args.bia_list) if v.amp == amp]
    
    # Make array of all bias images for given amp
    array_images = np.array([bia.image for bia in args.bia_list[sel]])
    
    # Measure the biweight midvariance (sigma) for a given pixel and take
    # the biweight average over all sigma to reduce the noise in the first 
    # measurement.
    S = biweight_location(biweight_midvariance(array_images,axis=(0,)))
    
    return S
    

def measure_gain(args, amp, rdnoise, flow=500, fhigh=35000, fnum=35):
    sel = [i for i,v in enumerate(args.ptc_list) if v.amp == amp]
    s_sel = sel[np.array([args.ptc_list[i].basename for i in sel]).argsort()]
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
        array_diff[i,:,:] = F1*m1/m2 - F2
    bins = np.logspace(np.log10(flow), np.log10(fhigh), fnum)
    gn = []
    for i in xrange(len(bins)-1):
        yloc, xloc = np.where((array_avg>bins[i]) * (array_avg<bins[i+1]))
        std = biweight_midvariance(array_diff[yloc,xloc])
        vr   = (std**2 - 2.*rdnoise**2) / 2.
        mn = biweight_location(array_avg)
        print("%s | Gain: %01.3f | RDNOISE: %01.3f | <ADU>: %0.1f | Pixels: %i" 
                %(amp, mn / vr, mn / vr * rdnoise, mn, len(xloc))) 
        gn.append(mn/vr)
    return biweight_location(gn)


def make_pixelflats(args):
    pass
    

def relative_throughput(args):
    pass


def main():
    # Read the arguments from the command line
    args = parse_args()
    
    # Get the bias jumps/structure for each amp
    (biasjump_left, biasjump_right, structure, 
     overscan, masterbias) = {}, {}, {}, {}, {}
    progress = ProgressBar(len(AMPS), 'Checking Biases %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        (biasjump_left[amp], biasjump_right[amp], 
         structure[amp], overscan[amp], 
         masterbias[amp]) = check_bias(args, amp)
        progress.current+=1
        progress()
    progress.done()
    
    # Get the dark jumps/structure and average counts
    darkcounts, masterdark = {}, {}
    progress = ProgressBar(len(AMPS), 'Checking Darks %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        darkcounts[amp], masterdark[amp] = check_darks(args, amp, 
                                                       masterbias[amp])
        progress.current+=1
        progress()
    progress.done()
    
    # Get the readnoise for each amp
    readnoise = {}
    progress = ProgressBar(len(AMPS), 'Measuring Readnoise %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        readnoise[amp] = measure_readnoise(args, amp)
        progress.current+=1
        progress()
    progress.done()
    
    # Get the readnoise for each amp
    gain = {}
    progress = ProgressBar(len(AMPS), 'Measuring Gain %s' %args.specid, 
                           fmt=ProgressBar.FULL)
    for amp in AMPS:
        gain[amp] = measure_gain(args, amp, readnoise[amp])
        progress.current+=1
        progress()
    progress.done()
        
if __name__ == '__main__':
    main()  