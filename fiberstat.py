# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 12:04:51 2017

@author: gregz
"""

import numpy as np
from utils import biweight_location, matrixCheby2D_7
import argparse as ap
import os.path as op
from astropy.io import fits
import sys
import logging
import matplotlib.pyplot as plt

def setup_logging():
    log = logging.getLogger('panacea')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        level = logging.INFO
           
        fmt = logging.Formatter(fmt)
        
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
        
        log = logging.getLogger('fiberstat')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log
    
def parse_args(argv=None):
    # Arguments to parse include ssp code, metallicity, isochrone choice, 
    #   whether this is real data or mock data
    parser = ap.ArgumentParser(description="MCSED",
                            formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-f","--filename", 
                        help='''Sky subtracted frame to be read, "S*.fits"''',
                        type=str, default=None)
                            
    parser.add_argument("-d","--dist", 
                        help='''Corresponding distortion file, ".dist"''',
                        type=str, default=None)

    parser.add_argument("-o","--outfile", 
                        help='''Name of the outfile''',
                        type=str, default=None)

    parser.add_argument("-i","--image_seg", 
                        help='''Segment of image to examine, if None, use all
                        Ex. "500,700,1200,1400" which is "xl,xh,yl,yh"''',
                        type=str, default=None)
                   
    args = parser.parse_args(args=argv)
    if args.image_seg is not None:
        args.image_seg = [int(i) for i in args.image_seg.split(',')]
    return args

def Fiberstat(image, D, outname, log, im_seg, fbins=12, fmax=6.):
    log.info('Creating fiber distance matrix')

    Y, X = np.indices(image.shape)
    Vxy = matrixCheby2D_7(D._scal_x(X.ravel()), 
                             D._scal_y(Y.ravel()))
    F = np.dot(Vxy, D.fiber_par_.data)
    F = F.reshape(image.shape)
    F0 = D.reference_f_.data
    offsetsbyfiber = np.abs(F[:,:,np.newaxis] - F0)
    fdist = np.min(offsetsbyfiber, axis=2)
    
    frange = np.linspace(0,fmax,fbins+1)
    
    if im_seg is not None:
        fdist = fdist[im_seg[2]:im_seg[3],im_seg[0]:im_seg[1]]
        image = image[im_seg[2]:im_seg[3],im_seg[0]:im_seg[1]]
    fdist = fdist.ravel()
    image = image.ravel()
    stats = np.zeros((fbins,))
    n = np.zeros((fbins,))
    
    log.info('Calculating average value as a function of fiber distance')
    
        
    for i in xrange(fbins):
        sel = np.where((fdist>=frange[i])  * (fdist<frange[i+1]))[0]
        stats[i] = biweight_location(image[sel])
        n[i] = len(sel)
    plt.figure(figsize=(6,5))
    x = frange[:-1]+np.diff(frange)/2.
    
    log.info('Plotting results')
    
    plt.plot(x , stats, color=[0./255.,175./255.,202./255.], lw=2.2)
    for xi, ni, si in zip(x, n, stats):
        plt.text(xi, si-0.04, '%0.1e' %ni, fontsize=8)
    plt.xlim([0, fmax])
    plt.xlabel('Fiber Distance', fontsize=14)
    plt.ylabel('Average Value', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(outname, dpi=150)
    plt.close()
    
def main():
    args = parse_args()
    
    log = setup_logging()
    
    try:
        from pyhetdex.cure.distortion import Distortion
        log.info('pyhetdex successfully loaded')
    except:
        log.error(''''No pyhetdex found.  Please install from: pip install --update --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex''')
        sys.exit(1)
    
    D = Distortion(args.dist)
    image = fits.open(args.filename)[0].data
    if args.outfile is None:
        args.outfile = 'fiberstat_%s.png' %(op.basename(args.filename)[:-5])
    Fiberstat(image, D, args.outfile, log, args.image_seg)
    
if __name__ == '__main__':
    main()