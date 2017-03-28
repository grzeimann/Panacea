# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 08:49:32 2017

@author: gregz
"""

import os.path as op
import logging
from astropy.io import fits
import numpy as np
from distutils.dir_util import mkpath

__all__ = ["Spectrograph"]

class Spectrograph:
    def __init__(self, path, specid, ifuslot, ifuid, basename, ifucen_fn=None,
                 verbose=True, N=None, D=None):
        self.path = path
        self.specid = specid
        self.ifuslot = ifuslot
        self.ifuid = ifuid
        self.basename = basename
        self.verbose = verbose
        self.setup_logging()
        self.ifucen_fn = ifucen_fn
        self.side_dict = {"L": ["LL","LU"], "R": ["RU","RL"]}
        self.N = N
        self.D = D

    def setup_logging(self):
        '''Set up a logger for shuffle with a name ``spectrograph``.
    
        Use a StreamHandler to write to stdout and set the level to DEBUG if
        verbose is set from the command line
        '''
        self.log = logging.getLogger('panacea')
        if not len(self.log.handlers):
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            if not self.verbose:
                level = logging.WARNING
            else:
                level = logging.INFO
           
            fmt = logging.Formatter(fmt)
        
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)
        
            self.log = logging.getLogger('panacea')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

    def write_to_fits(self, hdu, outname):
        '''
        Writing fits file to outname
        '''
        mkpath(op.dirname(outname))
        try:
            hdu.writeto(outname, overwrite=True)
        except TypeError:
            hdu.writeto(outname, clobber=True)
            
    def write_spectrograph_image(self, side, ext, prefix):
        if side is not "L" and side is not "R":
            self.log.error('The specified side needs to "L" or "R"')
            return None
        image = []
        outname = op.join(self.path,'%s%s_%s_sci_%s.fits' 
                                    %(prefix, self.basename.split('_',[]), 
                                      self.ifuslot, side))
                                  
        self.log.info('Makeing spectrograph image %s' %op.basename(outname))
        for i, amp in enumerate(self.side_dict[side]):
            fn = op.join(self.path, 'multi_%s_%s_%s_%s.fits' %(self.specid,
                                                               self.ifuslot,
                                                               self.ifuid,
                                                               amp))
            try:
                F = fits.open(fn)
                image.append(F[ext].data)
            except IOError:
                self.log.error('Failed to open %s' %fn)
                if self.N is None:
                    self.log.error('Need to know the size of the image to make'
                                   ' it all zeros.')
                    return None
                else:
                    image.append(np.zeros((self.N,self.D)))
        new = np.zeros((self.N*2, self.D))
        new[:self.N,:] = image[0]
        new[self.N:,:] = image[1]
        hdu = fits.PrimaryHDU(np.array(new, dtype='float32'))
        hdu.header.remove('BIASSEC')
        hdu.header.remove('TRIMSEC')
        hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,self.D,1,2*self.N)
        self.write_to_fits(hdu, outname)
                    
                
                
            
            