# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 19:14:32 2017

@author: gregz
"""

from astropy.io import fits
import glob
import os.path as op
from distutils.dir_util import mkpath
import argparse as ap
import numpy as np
import re


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
    
    description = " Prep VIRUS-W for Panacea Reductions "
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("--folder", nargs='?', type=str,
                        help='''Folder with raw virus-w data''')
                        
    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory
                        Default: \"virusw_raw"''', 
                        default="virusw_raw")
                        
    parser.add_argument("--overscan_pixel_length", nargs='?', type=int,
                        help='''number of pixels in overscan''',
                        default=50)                        
                        
    args = parser.parse_args(args=argv)

    return args
    
def check_criteria(header):
    """Check file if we want the file
    
    Parameters
    ----------
    header : fits header object
        from astropy.io.fits.open function
        
    Returns
    -------
    Bool
        True or False
    """
        
    image_list = ['sky', 'object']
    if header['imagetyp'] in image_list:
        return True, 'sci'
        
    object_list = ['twighlight']
    if header['object'] in object_list:
        return True, 'twi'
        
def build_fits(image, args, half, imtype, date, exptime):
    """Build fits object with proper header for Panacea
    
    Parameters
    ----------
    image : numpy array
        fits image
    args : NameSpace
        parsed args
    half : str
        'L' for bottom and 'U' for top
    imtype : str
        ex: 'sci' or 'twi' 
    date : str
        ex: '2017-1-03'
    exptime : float
        ex: 360. 
        
    Returns
    -------
    PrimaryHDU Object
        fits object to be written later
    """
    
    image = np.rot90(image)    
    hdu = fits.PrimaryHDU(image)
    a,b = image.shape
    if imtype == 'L':
        x1,x2,y1,y2 = (1, b, a + 1 - args.overscan_pixel_length, a)
        hdu.header['BIASSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
        x1,x2,y1,y2 = (1, b, 1, a - args.overscan_pixel_length)
        hdu.header['TRIMSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
    else:
        x1,x2,y1,y2 = (1, b, 1, args.overscan_pixel_length)
        hdu.header['BIASSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
        x1,x2,y1,y2 = (1, b, args.overscan_pixel_length + 1, a)
        hdu.header['TRIMSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
    hdu.header['GAIN'] = 1.0
    hdu.header['RDNOISE'] = 1.0
    hdu.header['CCDPOS'] = 'L'
    hdu.header['CCDHALF'] = half
    hdu.header['IMAGETYP'] = imtype    
    hdu.header['SPECID'] = 000
    hdu.header['IFUSLOT'] = 000
    hdu.header['IFUID'] = 000
    hdu.header['DATE-OBS'] = date
    hdu.header['EXPTIME'] =exptime
    
    return hdu
    
def write_to_fits(hdu, outname):
    """Write fits out

    Parameters
    ----------
    hdu : PrimaryHDU Object
        fits obect to write out
        
    outname : str
        file name to write out
    """
    
    try:
        hdu.writeto(outname, overwrite=True)
    except TypeError:
        hdu.writeto(outname, clobber=True)
    
def main():
    """
    
    Main Loop
    
    """
    args = parse_args()
    file_list = glob.glob(op.join(args.folder, '*.fits'))
    obsid = 0
    exp_num = 1
    unique_objects = []
    for file_name in file_list:
        F = fits.open(file_name)
        check, imtype = check_criteria(F[0].header)
        if check:
            a,b  = F[0].data.shape
            date = F[0].header['DATE']
            datefolder = ''.join(re.split('-',F[0].header['DATE-OBS']))
            object_name = F[0].header['OBJECT'].split(' ')[0]
            if object_name not in unique_objects:
                unique_objects.append(object_name)
                obsid +=1
                exp_num = 1
            exptime = F[0].header['EXPTIME']
            datetime = ''.join(F[0].header['TIME'].split(' ')[0].split(':'))[:-1]
            half = ['L', 'U']
            for h in half:
                F1 = build_fits(F[0].data[:,:int(b/2)], args, h, imtype, 
                            date, exptime)
                path = op.join(args.output,datefolder, 'virusw', 
                               'virusw%07d' %obsid, 'exp%02d' %exp_num, 
                               'virusw')
                mkpath(path)
                outname = op.join(path, '%sT%s_000L%s_%s.fits' 
                                  %(datefolder, datetime, h, imtype))
                write_to_fits(F1, outname)
            
            exp_num += 1
                                