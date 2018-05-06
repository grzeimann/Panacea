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
from utils import biweight_location, matrixCheby2D_7

__all__ = ["Spectrograph"]

class Spectrograph:
    def __init__(self, path, specid, ifuslot, ifuid, basename, ifucen=None,
                 verbose=True, N=None, D=None, nfib=None, wavelim=None,
                 disp=None, scale=None, collapselim=None, side_dict=None,
                 sides=None, header=None):
        self.path = path
        self.specid = specid
        self.ifuslot = ifuslot
        self.ifuid = ifuid
        self.basename = basename
        self.verbose = verbose
        self.setup_logging()
        self.ifucen = ifucen
        self.side_dict = side_dict
        self.sides = sides
        self.N = N
        self.D = D
        self.nfib = nfib
        self.wavelim = wavelim
        self.collapselim = collapselim
        self.disp = disp
        self.scale = scale
        self.fiberextract = {"L": None, "R": None}
        self.seeing = scale * 1.5
        self.header = header
        
    def setup_logging(self):
        '''Set up a logger for shuffle with a name ``panacea``.
    
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
    
    def check_side(self, side):
        if side is not "L" and side is not "R":
            self.log.error('The specified side needs to "L" or "R"')
            return False
        else:
            return True
            
    def build_outname(self, side, prefix):
        return op.join(self.path,'%s%s_%s_sci_%s.fits' 
                                 %(prefix, self.basename.split('_')[0], 
                                   self.ifuslot, side))

    def build_cubename(self, prefix):
        return op.join(self.path,'%s%s_%s_sci.fits' 
                                 %(prefix, self.basename.split('_')[0], 
                                   self.ifuslot))

    def build_inputname(self, amp):
        return op.join(self.path, 'multi_%s_%s_%s_%s.fits' %(self.specid,
                                                             self.ifuslot,
                                                             self.ifuid,
                                                             amp))
                                                             
    def load_file(self, amp):
        fn = self.build_inputname(amp)
        try:
            F = fits.open(fn)
            return F
        except IOError:
            self.log.error(
            'Failed to open %s' %fn)
            return None
            
    def write_spectrograph_image(self, side, ext, prefix):
        if not self.check_side(side):
            return None
        outname = self.build_outname(side, prefix)
        self.log.info('Making spectrograph image for %s' %op.basename(outname))
        image = []
        for i, amp in enumerate(self.side_dict[side]):
            F = self.load_file(amp)
            if F is not None:
                image.append(F[ext].data)
            else:
                image.append(np.zeros((self.N, self.D)))
        new = np.zeros((self.N*2, self.D))
        new[:self.N,:] = image[0]
        new[self.N:,:] = image[1]
        hdu = fits.PrimaryHDU(np.array(new, dtype='float32'), 
                              header=self.header)
        hdu.header.remove('BIASSEC')
        hdu.header.remove('TRIMSEC')
        hdu.header.remove('AMPSEC')
        hdu.header.remove('DETSEC')
        hdu.header.remove('CCDSEC')
        hdu.header['DATASEC']='[%i:%i,%i:%i]'%(1,self.D,1,2*self.N)
        hdu.header['CCDPOS']=side
        hdu.header.remove('CCDHALF')
        hdu.header.remove('AMPLIFIE')
        self.write_to_fits(hdu, outname)


    def build_FE(self, side, ext):
        image = []
        wave = []
        size=0
        for i, amp in enumerate(self.side_dict[side]):
            F = self.load_file(amp)
            if F is not None:
                image.append(F[ext].data)
                size += F[ext].data.shape[0]
                wave.append(F['wavelength'].data)
            else:
                image.append(np.zeros((112,self.D)))
                size += 112
                wave.append(np.zeros((112,self.D)))
        new = np.zeros((size, self.D))
        new[:image[1].shape[0],:] = image[1][::-1,:]
        new[image[1].shape[0]:,:] = image[0][::-1,:]
        newwave = np.zeros((size, self.D))
        newwave[:wave[1].shape[0],:] = wave[1][::-1,:]
        newwave[wave[1].shape[0]:,:] = wave[0][::-1,:]
        wv = np.arange(self.wavelim[0], self.wavelim[1]+self.disp, self.disp)
        newspec = np.zeros((size, len(wv)))
        for i in xrange(size):
            diff_wave = np.diff(newwave[i,:])
            wi = newwave[i,:-1]
            df = np.interp(wv, wi, diff_wave, left=0.0, right=0.0)
            fl = np.interp(wv, wi, new[i,:-1], left=0.0, right=0.0)
            newspec[i,:] = np.where(df!=0, fl/df*self.disp, 0.0)
        self.fiberextract[side] = newspec

    def write_fiberextract(self, side, ext, prefix):
        if not self.check_side(side):
            return None
        outname = self.build_outname(side, prefix)
        self.log.info('Making fiberextract image for  %s' %op.basename(outname))
        self.build_FE(side, ext)
        hdu = fits.PrimaryHDU(np.array(self.fiberextract[side], 
                                       dtype='float32'), header=self.header)
        hdu.header.remove('BIASSEC')
        hdu.header.remove('TRIMSEC')
        hdu.header.remove('AMPSEC')
        hdu.header.remove('DETSEC')
        hdu.header.remove('CCDSEC')
        hdu.header['CRVAL1'] = self.wavelim[0]
        hdu.header['CDELT1'] = self.disp
        hdu.header['CD1_1'] = self.disp
        hdu.header['CRPIX1'] = 1
        hdu.header['CCDPOS']=side
        hdu.header['DATASEC']='[%i:%i,%i:%i]'%(1,self.fiberextract[side].shape[0],
                                               1,self.fiberextract[side].shape[1])
        hdu.header.remove('CCDHALF')
        hdu.header.remove('AMPLIFIE')
        self.write_to_fits(hdu, outname)  

    def write_cube(self, ext, prefix, side=None):
        if side is not None:
            outname = self.build_outname(side, prefix[0])
            outname2 = self.build_outname(side, prefix[1])
            side_list = [side]
        else:
            outname = self.build_cubename(prefix[0])
            outname2 = self.build_cubename(prefix[1])
            side_list = self.sides
        self.log.info('Making cube image for %s' %op.basename(outname))
        data = []
        for side in side_list:
            if self.fiberextract[side] is None:
                self.build_FE(side, ext)
            data.append(self.fiberextract[side])
        data = np.vstack(data)
        a,b = data.shape
        x = np.arange(self.ifucen[:,1].min()-self.scale, 
                      self.ifucen[:,1].max()+self.scale, self.scale)
        y = np.arange(self.ifucen[:,2].min()-self.scale, 
                      self.ifucen[:,2].max()+self.scale, self.scale)
        xgrid, ygrid = np.meshgrid(x, y)
        zgrid = np.zeros((b,)+xgrid.shape)
        d = np.zeros((len(self.ifucen[:,1]),)+xgrid.shape)
        w = np.zeros((len(self.ifucen[:,1]),)+xgrid.shape)
        for i in xrange(len(x)):
            for j in xrange(len(y)):
                d[:,j,i]= np.sqrt((self.ifucen[:,1] - xgrid[j,i])**2 + 
                            (self.ifucen[:,2] - ygrid[j,i])**2)
                w[:,j,i] = np.exp(-1./2.*(d[:,j,i]/self.seeing)**2)    
        ws = w.sum(axis=0)
        for k in xrange(b):
            zgrid[k,:,:] = (data[:,k][:,np.newaxis,np.newaxis]*w).sum(axis=0)/ws
        
                    
        hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
        hdu.header['CRPIX1'] = x[0]
        hdu.header['CRPIX2'] = y[0]
        hdu.header['CDELT3'] = self.disp
        hdu.header['CRVAL3'] = self.wavelim[0]
        hdu.header['CRPIX3'] = 1
        self.write_to_fits(hdu, outname)
        self.log.info('Making collapsed cube image for %s' %op.basename(outname2))
        wv = np.arange(self.wavelim[0], self.wavelim[1]+self.disp, self.disp)
        sellow = np.searchsorted(wv, self.collapselim[0], side='left')
        selhigh = np.searchsorted(wv, self.collapselim[1], side='right')
        zimage = biweight_location(zgrid[sellow:selhigh,:,:],axis=(0,))
        hdu = fits.PrimaryHDU(np.array(zimage, dtype='float32'))
        hdu.header['CRPIX1'] = x[0]
        hdu.header['CRPIX2'] = y[0]
        self.write_to_fits(hdu, outname2)
        
        
    def write_new_distfile(self, D, side):
        outname = op.join(self.path,'mastertrace_%s_%s.dist' 
                                 %(self.specid, side))
        self.log.info('Making distortion frame for %s' %op.basename(outname))
        trace = []
        wave = []
        size=0
        for i, amp in enumerate(self.side_dict[side]):
            F = self.load_file(amp)
            if F is not None:
                trace.append(F['trace'].data)
                size += F['trace'].data.shape[0]
                wave.append(F['wavelength'].data)
            else:
                trace.append(np.zeros((self.nfib,self.D)))
                size += self.nfib
                wave.append(np.zeros((self.nfib,self.D)))
        col = int(self.D / 2)
        intv = [1, 1+self.N]
        ypos = np.zeros((size, self.D))
        ypos[:trace[1].shape[0],:] = trace[1][::-1,:]+intv[1]
        ypos[trace[1].shape[0]:,:] = trace[0][::-1,:]+intv[0]
        xpos = np.ones((size,1)) * np.arange(self.D)
        f1 = ypos[:,col]
        f0 = np.sort(f1)[::-1]
        fpos = np.zeros((ypos.shape))
        for k in np.arange(size):
            fpos[k,:] = f1[k]
        wpos = np.zeros((size, self.D))
        wpos[:wave[1].shape[0],:] = wave[1][::-1,:]
        wpos[wave[1].shape[0]:,:] = wave[0][::-1,:]      
        
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
        D.writeto(outname)

