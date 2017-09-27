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
from utils import biweight_location


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

    parser.add_argument("-hi","--highres", help='''High Res Mode?''',
                        action="count", default=0)
                        
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
    try:
        if header['IMAGETYP'] in image_list:
            return True, 'sci'
    except:
        return False, ''
        
    object_list = ['twilight']
    try:
        if header['OBJECT'] in object_list:
            return True, 'twi'
    except:
        return False, ''
 
    object_list = ['zero']
    try:
        if header['IMAGETYP'] in object_list:
            return True, 'zro'
    except:
        return False, ''
   
    return False, ''
        
def build_fits(image, args, half, side, imtype, date, exptime, object_name,
               file_name):
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
    image = np.fliplr(image)    
    a,b = image.shape        
    hdu = fits.PrimaryHDU(image)
    x1,x2,y1,y2 = (1, 44, 1, a)
    hdu.header['BIASSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
    x1,x2,y1,y2 = (45, b, 1, a)
    hdu.header['TRIMSEC'] = '[%i:%i,%i:%i]' %(x1,x2,y1,y2)
    hdu.header['GAIN'] = 0.62
    hdu.header['RDNOISE'] = 1.0
    hdu.header['CCDPOS'] = side
    hdu.header['CCDHALF'] = half
    hdu.header['IMAGETYP'] = imtype    
    hdu.header['SPECID'] = 0
    hdu.header['IFUSLOT'] = 0
    hdu.header['IFUID'] = '000'
    hdu.header['DATE-OBS'] = date
    hdu.header['EXPTIME'] = exptime
    hdu.header['OBJECT'] = object_name
    hdu.header['FILENAME'] = file_name
    
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
    
    unique_date = {}
    
    file_list = sorted(glob.glob(op.join(args.folder, '*.fits')))
    for file_name in file_list:
        try:
            F = fits.open(file_name)    
            # Checking if it is a sci or twi, otherwise move on
            check, imtype = check_criteria(F[0].header)
        except:
            check = False
            imtype = ''
            print('Error reading in %s' %file_name)
        if check:
            # Bookkeeping for date/obs/expn
            datefolder = ''.join(F[0].header['DATE-OBS'].split('T')[0].split('-'))
            object_name = F[0].header['OBJECT'].split(' ')[0]
            if datefolder not in unique_date:
                unique_date[datefolder] = {}
                unique_date[datefolder][object_name] = [1,1]
            else:
                if object_name not in unique_date[datefolder]:
                    o_cnt = len(unique_date[datefolder]) + 1
                    unique_date[datefolder][object_name] = [o_cnt,1]
                else:
                    unique_date[datefolder][object_name][1] += 1
                    
            obsid = unique_date[datefolder][object_name][0]
            exp_num = unique_date[datefolder][object_name][1]
            print(file_name,datefolder, object_name, obsid, exp_num)

            a,b  = F[0].data.shape
            date = F[0].header['DATE']
            exptime = F[0].header['EXPTIME']
            datetime = ''.join(F[0].header['TIME'].split(' ')[0].split(':'))[:-1]
            if args.highres:
                side = 'R'
                half = 'U'
            else:
                side = 'L'
                half = 'L'
            data = np.hstack([F[0].data[:,:1024]
                              -biweight_location(F[0].data[-43:,:1024]),
                              F[0].data[:,1124:]
                              -biweight_location(F[0].data[-43:,1124:])])
            F1 = build_fits(data, args, half, side, imtype, 
                            date, exptime, object_name, op.basename(file_name))
            path = op.join(args.output, datefolder, 'virusw', 
                           'virusw%07d' %obsid, 'exp%02d' %exp_num, 
                           'virusw')
            mkpath(path)
            outname = op.join(path, '%sT%s_000%s%s_%s.fits' 
                              %(datefolder, datetime, side, 'U', imtype))
            write_to_fits(F1, outname)
            
                
if __name__ == '__main__':
    main()    
                                            