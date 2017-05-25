# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:55:42 2017

@author: gregz
"""

import argparse as ap
import glob
import os.path as op
from astropy.io import fits
from pyhetdex.het.fplane import FPlane
from pyhetdex.coordinates.tangent_projection import TangentPlane as TP
import numpy as np
from utils import biweight_location, biweight_filter
from scipy.signal import medfilt
import logging
import sys


# common parts of the argument parser
astro_parent = ap.ArgumentParser(add_help=False)
_group_astro = astro_parent.add_mutually_exclusive_group()
_group_astro.add_argument("--ifuslots", nargs='?', type=str, 
                            help='''IFUSLOT for processing. 
                            Ex: "075,076".''', default = None)
_group_astro.add_argument("--ifuslot_file", nargs='?', type=str, 
                            help='''IFUSLOT file. 
                            Ex: "good_ifuslots.txt".''', default = None)
_group_astro.add_argument("--folder", nargs='?', type=str, 
                            help='''Reduction Folder
                            Ex: "reductions".''', default = "reductions")
_group_astro.add_argument("--instr", nargs='?', type=str, 
                            help='''Instrument to process. 
                            Default: "virus"
                            Ex: "camra" for lab data,
                                "lrs2" for lrs2.''', default = "virus")
_group_astro.add_argument("-sd","--scidir_date", nargs='?', type=str,
                            help='''Science Directory Date.    
                            Ex: \"20160412\"''', default=None)
_group_astro.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                            help='''Science Directory ObsID.   
                            Ex: \"3\" or \"102\"''', default=None)                            
_group_astro.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                            help='''Science Directory exposure number.
                            Ex: \"1\" or \"05\"''', default=None)              
_group_astro.add_argument("--fplane", nargs='?', type=str, 
                            help='''Fplane file
                            Ex: "fplane.txt".''', default = 'fplane.txt')

def setup_logging():
    '''Set up a logger for analysis with a name ``shot``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('shot')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        level = logging.INFO
       
        fmt = logging.Formatter(fmt)
    
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
    
        log = logging.getLogger('shot')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log

def astrometry(opts, ifuslot, pos, ra, dec, pa, sys_rot=1.6):
    # Read in the Fplane File
    fplane = FPlane(opts.fplane)
    
    # Set up the Tangent Plane projection
    tp = TP(ra, dec, pa+sys_rot, 1., 1.)

    # Get IFU positions for given ifuslot
    ifu = fplane.by_ifuslot(ifuslot)

    # remember to flip x,y 
    xfp = pos[:,0] + ifu.y
    yfp = pos[:,1] + ifu.x
    
    ra_dec = np.zeros((len(xfp),2))
    ra_dec[:,0], ra_dec[:,1] = tp.xy2raDec(xfp, yfp)
    return ra_dec
 
def read_in_multifiles(args, ifuslot, extensions=None):
    '''
    Read in the multi_* files from panacea reduction
    
    Return arrays
    '''
    log = setup_logging()
    if extensions is None:
        log.error('Extensions were not provided for input')
        sys.exit(1)   
        
    # Reading in the scidir_date(obsid, expnum) provided
    obs = 'sci'                         
    labels = ['dir_date', 'dir_obsid', 'dir_expnum']

    for label in labels[:2]:
        if getattr(args, obs+label) is None:
            log.error('%s%s was not provided' %(obs, label))
            sys.exit(1) 
        if not isinstance(getattr(args, obs+label), list):
            setattr(args, obs+label, 
                    getattr(args, obs+label).replace(" ", "").split(','))
    if getattr(args, obs+labels[2]) is not None \
       and not isinstance(getattr(args, obs+label[2]), list):
        setattr(args, obs+labels[2], 
                getattr(args, obs+labels[2]).replace(" ", "").split(','))
    tuple_outer = ()
    for i, ext in enumerate(extensions):
        tuple_outer.append([])
    for date in getattr(args, obs+labels[0]):
        for obsid in getattr(args, obs+labels[1]):
            if getattr(args, obs+labels[2]) is not None:  
                exposures = getattr(args, obs+labels[2])
            else:
                folder = op.join(date, args.instr,
                                 "{:s}{:07d}".format(args.instr, 
                                                     int(obsid)))
                exposures = [op.basename(exp)[-2:] for exp in 
                                   glob.glob(op.join(args.folder, folder,'*'))] 
            for expnum in exposures:
                tuple_inner = ()
                for i, ext in enumerate(extensions):
                    tuple_inner.append([])
                folder = op.join(date, 
                                 args.instr, 
                                 "{:s}{:07d}".format(args.instr,
                                 int(obsid)), 
                                 "exp{:02d}".format(int(expnum)),
                                 args.instr)
                files = sorted(glob.glob(op.join(args.folder, folder, 
                                                 'multi_*_%s*' %ifuslot)))
                for fn in files:
                    F = fits.open(fn)
                    for i, ext in enumerate(extensions):
                        tuple_inner[i].append(F[ext].data)
                for i, ext in enumerate(extensions):
                    tuple_outer[i].append(np.vstack(tuple_inner[i]))
    return tuple_outer

def fiber_to_fiber_shot(args=None):
    '''
    Parameters
    ----------
    args : list of strings, optional
        command line
    '''
    parser = ap.ArgumentParser(description="""Calculate the fiber
                                     to fiber across a shot.""",
                                     parents=[astro_parent, ])
                                     
    args = parser.parse_args(args)
    
    for ifuslot in args.ifuslot:
        wave_array, twi_array = read_in_multifiles(args, ifuslot, 
                                                   extensions=['twi_spectrum', 
                                                               'wavelength'])
    std_wave = np.arange(wvl,wvh+wvstep,wvstep)
    a,b = twi_array.shape
    large_array = []
    for fib in np.arange(a):
        diff = np.diff(wave_array[fib,:])
        diff_array = np.hstack([diff,diff[-1]])
        y = np.interp(std_wave, wave_array[fib,:],
                      twi_array[fib,:], left=-999., right=-999.)
        large_array.append(y/diff_array) 
    average = np.vstack(large_array)
    average = np.ma.array(average, mask=(average==-999.))
    average = biweight_location(average,axis=(0,))
    fiber_to_fiber = []
    for fib in np.arange(a):
        y = medfilt(large_array[fib]/average, scale)
        fiber_to_fiber.append(np.interp(wave_array[fib,:],
                                        std_wave, y))
    pass
   
def collapse_fibers():
    '''
    spectrum, fiber_to_fiber (all), ifupos, grid schema, wave1, wave2
    '''
    # collapse the spectrum/fiber_to_fiber (flag and remove fibers that are zero
    #     fiber to fiber using an biweight over wave1, wave2
    # gaussian interpolate and conserve flux for each grid point (see spectrograph.py)
    # return image
    
    pass

def get_wcs_collapsed_image():
    '''
    Not sure about this yet, but certainly needed
    '''
    pass
    
def subtract_image_background():
    '''
    image, mask
    '''
    # Follow https://photutils.readthedocs.io/en/stable/photutils/background.html
    # Need to time this step quite a bit
    # iterate with detections?
    pass

def detect_in_image():
    pass

def match_to_database():
    '''
    RA, Decdetection catalogs and vizier (astroquery) match?  vo conesearch?    
    '''
    pass



def throughput(): 
    pass

def find_positional_offset():
    pass

def recall_astrometry_using_offset():
    pass
    
def visualize_sig_map():
    pass

