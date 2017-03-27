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
import time
import re

import numpy as np
import matplotlib.pyplot as plt
import os.path as op
from astropy.io import fits
from pyhetdex.cure.distortion import Distortion
from pyhetdex.het.ifu_centers import IFUCenter

from args import parse_args
from amplifier import Amplifier
from fiber_utils import get_model_image
from utils import matrixCheby2D_7, biweight_filter, biweight_midvariance
from utils import biweight_location
import config

import glob
                      

def imstat(image1, image2, Fiber1, Fiber2, outname, fbins=50, fmax=8.):
    a,b = image1.shape
    images = [image1,image2]
    Fibers = [Fiber1,Fiber2]
    totstat = np.zeros((2*a,b))
    totdist = np.zeros((2*a,b))
    for i,image in enumerate(images):
        dist = np.zeros((a,b))
        trace_array = np.array([fiber.trace for fiber in Fibers[i]])
        for y in xrange(a):
            for x in xrange(b):
                dist[y,x] = np.min(np.abs(trace_array[:,x] - y))
        if i==0:
            totdist[:a,:] = dist
            totstat[:a,:] = image
        else:
            totdist[a:,:] = dist
            totstat[a:,:] = image
    frange = np.linspace(0,fmax,fbins+1)
    totdist = totdist.ravel()
    totstat = totstat.ravel()
    stats = np.zeros((fbins,))
    for i in xrange(fbins):
        sel = np.where((totdist>=frange[i])  * (totdist<frange[i+1]))[0]
        stats[i] = biweight_location(totstat[sel])
    plt.figure(figsize=(6,5))
    plt.plot(frange[:-1]+np.diff(frange)/2., stats, color=[1.0, 0.2, 0.2], 
             lw=3)
    plt.xlim([0, fmax])
    plt.xlabel('Fiber Distance', fontsize=14)
    plt.ylabel('Average Value', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(outname, dpi=150)
    plt.close()
    
                      
def recreate_fiberextract(instr1, instr2, wavelim, disp):
    col = int(instr1.D / 2)
    intv = [1, 1+instr1.D]
    ypos = np.array([fiber.trace+intv[v] 
                     for v,instr in enumerate([instr1, instr2]) 
                     for fiber in instr.fibers])
    allspec = np.array([fiber.spectrum / fiber.fiber_to_fiber * (1-fiber.dead)
                        for v,instr in enumerate([instr1, instr2]) 
                        for fiber in instr.fibers])                                       
    allskys = np.array([(fiber.spectrum-fiber.sky_spectrum) 
                            / fiber.fiber_to_fiber * (1-fiber.dead)
                        for v,instr in enumerate([instr1, instr2]) 
                        for fiber in instr.fibers])
    allwave = np.array([fiber.wavelength 
                        for v,instr in enumerate([instr1, instr2]) 
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
    ypos = np.array([fiber.trace+intv[v] 
                     for v,instr in enumerate([instr1,instr2]) 
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
    
def make_cube_file(args, filename, ifucen, scale, side):
    if args.instr.lower() == "lrs2":
        outname = op.join(op.dirname(filename),'Cu' + op.basename(filename))
        outname2 = op.join(op.dirname(filename),'Co' + op.basename(filename))
        print("Making Cube image for %s" %op.basename(outname))
        try:
            F = fits.open(filename)
        except IOError:
            print("Could not open %s" %filename)
            return None
        data = F[0].data
        a,b = data.shape
        x = np.arange(ifucen[:,1].min()-scale, ifucen[:,1].max()+scale, scale)
        y = np.arange(ifucen[:,2].min()-scale, ifucen[:,2].max()+scale, scale)
        xgrid, ygrid = np.meshgrid(x, y)
        zgrid = np.zeros((b,)+xgrid.shape)
        for k in xrange(b):
            for i in xrange(len(x)):
                for j in xrange(len(y)):
                    d = np.sqrt((ifucen[:,1] - xgrid[j,i])**2 + 
                                (ifucen[:,2] - ygrid[j,i])**2)
                    w = np.exp(-1./2.*(d/1.)**2)
                    sel = w > 1e-3
                    zgrid[k,j,i] = np.sum(data[sel,k]*w[sel])/np.sum(w[sel])
        hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
        
        zcol = biweight_location(zgrid[int(b/3):int(2*b/3),:,:],axis=(0,))
        hdu.header['CDELT3'] = F[0].header['CDELT1']
        hdu.header['CRVAL3'] = F[0].header['CRVAL1']
        hdu.header['CRPIX3'] = F[0].header['CRPIX1']
        hdu.writeto(outname, overwrite=True)
        hdu = fits.PrimaryHDU(np.array(zcol, dtype='float32'))
        hdu.writeto(outname2, overwrite=True)
    if args.instr.lower() == "virus":
        if side == "R":
            file2 =filename
            file1 = filename[:-6] + "L.fits"
        else:
            return None
        if not op.exists(file2):
            print("Could not open %s" %file2)
            return None
        if not op.exists(file1):
            print("Could not open %s" %file1)
            return None
        outname = op.join(op.dirname(filename),'Cu' 
                                          + op.basename(filename)[:-7]+'.fits')
        outname2 = op.join(op.dirname(filename),'Co' 
                                          + op.basename(filename)[:-7]+'.fits')
        print("Making Cube image for %s" %op.basename(outname))
        F2 = fits.open(file2)
        F1 = fits.open(file1)
        data1 = F1[0].data
        data2 = F2[0].data
        a,b = data1.shape
        data = np.vstack([data1,data2])
        if len(data[:,0]) != len(ifucen[:,1]):
            print("Length of IFUcen file not the same as Fe. Skipping Cube")
            return None
        x = np.arange(ifucen[:,1].min()-scale, ifucen[:,1].max()+scale, scale)
        y = np.arange(ifucen[:,2].min()-scale, ifucen[:,2].max()+scale, scale)
        xgrid, ygrid = np.meshgrid(x, y)
        zgrid = np.zeros((b,)+xgrid.shape)
        for k in xrange(b):
            for i in xrange(len(x)):
                for j in xrange(len(y)):
                    d = np.sqrt((ifucen[:,1] - xgrid[j,i])**2 + 
                                (ifucen[:,2] - ygrid[j,i])**2)
                    w = np.exp(-1./2.*(d/(2./2.35*2.))**2)
                    sel = w > 1e-3
                    zgrid[k,j,i] = np.sum(data[sel,k]*w[sel])/np.sum(w[sel])
        hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
        zcol = biweight_location(zgrid[int(b/3):int(2*b/3),:,:],axis=(0,))
        hdu.header['CDELT3'] = F1[0].header['CDELT1']
        hdu.header['CRVAL3'] = F1[0].header['CRVAL1']
        hdu.header['CRPIX3'] = F1[0].header['CRPIX1']
        hdu.writeto(outname, overwrite=True)
        hdu = fits.PrimaryHDU(np.array(zcol, dtype='float32'))
        hdu.writeto(outname2, overwrite=True)
    
def make_error_frame(image1, image2, mask1, mask2, header, outname):
    print("Making error image for %s" %op.basename(outname))
    a,b = image1.shape
    new = np.zeros((a*2,b))
    mas = np.zeros((a*2,b))

    new[:a,:] = image1
    new[a:,:] = image2
    mas[:a,:] = mask1
    mas[a:,:] = mask2
    err = np.zeros(new.shape)
    for i in xrange(2*a):
        err[i,:] = np.where(mas[i,:]<0, mas[i,:], 
                            biweight_filter(new[i,:], 31, 
                                            func=biweight_midvariance))
    hdu = fits.PrimaryHDU(np.array(err, dtype='float32'), header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,2*a)
    outname = op.join(op.dirname(outname), 'e.' + op.basename(outname))
    hdu.writeto(outname, overwrite=True)   
    
def make_fiber_image(Fe, header, outname, args, amp):
    print("Making Fiberextract image for %s" %op.basename(outname))
    a,b = Fe.shape
    hdu = fits.PrimaryHDU(np.array(Fe, dtype='float32'), header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,a)
    hdu.header['CRVAL1'] = args.wvl_dict[amp][0]
    hdu.header['CDELT1'] = args.disp[amp]
    hdu.header['CD1_1'] = args.disp[amp]
    hdu.header['CRPIX1'] = 1
    hdu.writeto(outname, overwrite=True)    

def make_fiber_error(Fe, header, outname, args, amp):
    print("Making Fiberextract image for %s" %op.basename(outname))
    a,b = Fe.shape
    err = np.zeros(Fe.shape)
    for i in xrange(a):
        err[i,:] = biweight_filter(Fe[i,:], 21, func=biweight_midvariance)
    hdu = fits.PrimaryHDU(np.array(err, dtype='float32'), header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,a)
    hdu.header['CRVAL1'] = args.wvl_dict[amp][0]
    hdu.header['CDELT1'] = args.disp[amp]
    hdu.header['CD1_1'] = args.disp[amp]
    hdu.header['CRPIX1'] = 1
    outname = op.join(op.dirname(outname), 'e.' + op.basename(outname))
    hdu.writeto(outname, overwrite=True)  
    
def make_spectrograph_image(image1, image2, header, outname):
    print("Making spectrograph image for %s" %op.basename(outname))
    a,b = image1.shape
    new = np.zeros((a*2,b))
    new[:a,:] = image1
    new[a:,:] = image2
    hdu = fits.PrimaryHDU(np.array(new, dtype='float32'), header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,2*a)
    hdu.writeto(outname, overwrite=True)    

def make_amplifier_image(image, header, outname):
    print("Making amplifier image for %s" %op.basename(outname))
    a,b = image.shape
    hdu = fits.PrimaryHDU(np.array(image, dtype='float32'), header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,b,1,a)
    hdu.writeto(outname, overwrite=True) 
            
def reduce_science(args):
    for spec in args.specid:
        spec_ind_sci = np.where(args.sci_df['Specid'] == spec)[0]
        for amp in config.Amps:
            amp_ind_sci = np.where(args.sci_df['Amp'] == amp)[0]
            sci_sel = np.intersect1d(spec_ind_sci, amp_ind_sci) 
            for ind in sci_sel:
                if args.instr == "virus":
                    if not args.use_trace_ref:
                        ifucen = np.loadtxt(op.join(args.configdir, 
                                                    'IFUcen_files', 
                                                    args.ifucen_fn[amp][0]
                                                    + args.sci_df['Ifuid'][ind] 
                                                    + '.txt'), 
                                                    usecols=[0,1,2,4], 
                                               skiprows=args.ifucen_fn[amp][1])
                        
                    else:
                        if args.sci_df['Ifuid'][ind] == '004':
                            ifucen = np.loadtxt(op.join(args.configdir,
                                                    'IFUcen_files',
                                                    'IFUcen_HETDEX_reverse_R.txt'),
                                                    usecols=[0,1,2,4],
                                               skiprows=args.ifucen_fn[amp][1])
                            ifucen[224:,:] = ifucen[-1:223:-1,:]
                        else:
                            ifucen = np.loadtxt(op.join(args.configdir,
                                                    'IFUcen_files',
                                                    'IFUcen_HETDEX.txt'),
                                                    usecols=[0,1,2,4],
                                               skiprows=args.ifucen_fn[amp][1])
                else:
                    ifucen = np.loadtxt(op.join(args.configdir, 'IFUcen_files', 
                                        args.ifucen_fn[amp][0]), 
                              usecols=[0,1,2], skiprows=args.ifucen_fn[amp][1])
                if args.debug:
                    print("Working on Sci for %s, %s" %(spec, amp)) 
                if args.check_if_twi_exists:
                    fn = op.join(args.twi_dir,'fiber_*_%s_%s_%s_%s.pkl' %(spec, 
                                                   args.sci_df['Ifuslot'][ind],
                                                     args.sci_df['Ifuid'][ind],
                                                                          amp))
                    calfiles = glob.glob(fn)
                    if not calfiles:
                        print("No cals found for %s,%s: %s"
                              %(spec, amp, args.sci_df['Files'][ind]))
                        print("If you want to produce cals include "
                              "--reduce_twi")
                  
                sci1 = Amplifier(args.sci_df['Files'][ind],
                                 args.sci_df['Output'][ind],
                                 calpath=args.twi_dir, skypath=args.sky_dir,
                                 debug=False, refit=False, 
                                 dark_mult=args.dark_mult[amp],
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 use_trace_ref=args.use_trace_ref,
                                 calculate_shift=args.adjust_trace,
                                 fiber_date=args.fiber_date,
                                 cont_smooth=args.cont_smooth)
                #sci1.load_fibers()
                #if sci1.fibers and not args.start_from_scratch:
                #    if sci1.fibers[0].spectrum is not None:
                #        sci1.prepare_image()
                #        sci1.sky_subtraction()
                #        sci1.clean_cosmics()
                #else:
                sci1.load_all_cal()
                if args.adjust_trace:
                    sci1.refit=True
                    sci1.get_trace()
                    sci1.refit=False
                sci1.fiberextract()
                if args.refit_fiber_to_fiber:
                    sci1.refit=True
                    sci1.get_fiber_to_fiber()
                    sci1.refit=False
                sci1.sky_subtraction()
                sci1.clean_cosmics()
                sci1.fiberextract()
                sci1.sky_subtraction()
                sci2 = Amplifier(args.sci_df['Files'][ind].replace(amp, 
                                                      config.Amp_dict[amp][0]),
                                 args.sci_df['Output'][ind],
                                 calpath=args.twi_dir, skypath=args.sky_dir, 
                                 debug=False, refit=False, 
                             dark_mult=args.dark_mult[config.Amp_dict[amp][0]],
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 use_trace_ref=args.use_trace_ref,
                                 calculate_shift=args.adjust_trace,
                                 fiber_date=args.fiber_date,
                                 cont_smooth=args.cont_smooth)
                #sci2.load_fibers()
                #if sci2.fibers and not args.start_from_scratch:
                #    if sci2.fibers[0].spectrum is not None:
                #        sci2.prepare_image()
                #        sci2.sky_subtraction()
                #        sci2.clean_cosmics()
                #else:
                sci2.load_all_cal()
                if args.adjust_trace:
                    sci2.refit=True
                    sci2.get_trace()
                    sci2.refit=False
                sci2.fiberextract()
                if args.refit_fiber_to_fiber:
                    sci2.refit=True
                    sci2.get_fiber_to_fiber()
                    sci2.refit=False
                sci2.sky_subtraction()
                sci2.clean_cosmics()
                sci2.fiberextract()
                sci2.sky_subtraction()
                outname = op.join(args.sci_df['Output'][ind],
                                  'S%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_spectrograph_image(sci1.clean_image, sci2.clean_image, 
                                        sci1.header, outname)
                make_spectrograph_image(sci1.error, sci2.error, 
                                        sci1.header, op.join(op.dirname(outname), 'ee.'+op.basename(outname)))                   
                make_error_frame(sci1.clean_image, sci2.clean_image, sci1.mask,
                                 sci2.mask, sci1.header, outname)
                outname = op.join(args.sci_df['Output'][ind],
                                  'cS%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_spectrograph_image(np.where(sci1.mask==0, 
                                                 sci1.clean_image, 0.0),
                                        np.where(sci2.mask==0, 
                                                 sci2.clean_image, 0.0),
                                        sci1.header, outname)
                make_error_frame(sci1.clean_image, sci2.clean_image, sci1.mask,
                                 sci2.mask, sci1.header, outname)
                outname = op.join(args.sci_df['Output'][ind],
                                  'CsS%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_spectrograph_image(sci1.continuum_sub, sci2.continuum_sub, 
                                        sci1.header, outname)
                make_error_frame(sci1.continuum_sub, sci2.continuum_sub, 
                                 sci1.mask, sci2.mask, sci1.header, outname)
                outname = op.join(args.sci_df['Output'][ind],
                                  'cCsS%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_spectrograph_image(np.where(sci1.mask==0, 
                                                 sci1.continuum_sub, 0.0),
                                        np.where(sci2.mask==0, 
                                                 sci2.continuum_sub, 0.0),
                                        sci1.header, outname)
                make_error_frame(sci1.continuum_sub, sci2.continuum_sub,
                                 sci1.mask, sci2.mask, sci1.header, outname)
                outname = op.join(args.sci_df['Output'][ind],
                                  'cCsS%s_%s_sci_%s_imstat.png' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                imstat(sci1.residual, sci2.residual, sci1.fibers,
                       sci2.fibers, outname)
                Fe, FeS = recreate_fiberextract(sci1, sci2, 
                                                wavelim=args.wvl_dict[amp], 
                                                disp=args.disp[amp])
                outname = op.join(args.sci_df['Output'][ind],
                                  'Fe%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_fiber_image(Fe, sci1.header, outname, args, amp)
                make_fiber_error(Fe, sci1.header, outname, args, amp)
                make_cube_file(args, outname, ifucen, args.cube_scale, 
                               config.Amp_dict[amp][1])
                outname = op.join(args.sci_df['Output'][ind],
                                  'FeS%s_%s_sci_%s.fits' %(
                          op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                   args.sci_df['Ifuslot'][ind], 
                                                      config.Amp_dict[amp][1]))
                make_fiber_image(FeS, sci1.header, outname, args, amp)
                make_fiber_error(FeS, sci1.header, outname, args, amp)
                make_cube_file(args, outname, ifucen, args.cube_scale, 
                               config.Amp_dict[amp][1])
                if args.save_sci_fibers:
                    sci1.save_fibers()
                    sci2.save_fibers()
                if args.save_sci_amplifier:
                    sci1.save()
                    sci2.save()
                if args.debug:
                    print("Finished working on Sci for %s, %s" %(spec, amp))                    
                

def reduce_twighlight(args):
    D = Distortion(op.join(args.configdir, 'DeformerDefaults', 
                                        'mastertrace_twi_027_L.dist'))   
    for spec in args.specid:
        spec_ind_twi = np.where(args.twi_df['Specid'] == spec)[0]
        for amp in config.Amps:
            amp_ind_twi = np.where(args.twi_df['Amp'] == amp)[0]
            twi_sel = np.intersect1d(spec_ind_twi, amp_ind_twi)
            for ind in twi_sel:
                if args.debug:
                    print("Working on Cal for %s, %s" %(spec, amp))                    
                twi1 = Amplifier(args.twi_df['Files'][ind],
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True, dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 init_lims=args.wvl_dict[amp], 
                                 check_fibermodel=True, check_wave=True,
                                 fsize=args.fsize, 
                                 fibmodel_nbins=args.fibmodel_bins,
                                 sigma=args.fibmodel_sig,
                                 power=args.fibmodel_pow,
                                 use_trace_ref=args.use_trace_ref,
                                 default_fib = args.default_fib,
                                 wave_nbins = args.wave_nbins)
                #twi1.load_fibers()
                twi1.get_fiber_to_fiber()
                twi1.sky_subtraction()
                twi2 = Amplifier(args.twi_df['Files'][ind].replace(amp, 
                                                      config.Amp_dict[amp][0]),
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True,  dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 init_lims=args.wvl_dict[amp], 
                                 check_fibermodel=True, check_wave=True,
                                 fsize=args.fsize, 
                                 fibmodel_nbins=args.fibmodel_bins,
                                 sigma=args.fibmodel_sig,
                                 power=args.fibmodel_pow,
                                 use_trace_ref=args.use_trace_ref,
                                 default_fib = args.default_fib,
                                 wave_nbins = args.wave_nbins)
                #twi2.load_fibers()
                twi2.get_fiber_to_fiber()
                twi2.sky_subtraction()
                image1 = get_model_image(twi1.image, twi1.fibers, 
                                         'fiber_to_fiber', debug=twi1.debug)
                image2 = get_model_image(twi2.image, twi2.fibers, 
                                         'fiber_to_fiber', debug=twi2.debug)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'mastertrace_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    config.Amp_dict[amp][1]))
                make_spectrograph_image(image1, image2, twi1.header, outname)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'mastertwi_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    config.Amp_dict[amp][1]))  
                make_spectrograph_image(twi1.image, twi2.image, 
                                        twi1.header, outname)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'normtwi_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    amp))  
                make_amplifier_image(np.where(
                                 np.isfinite(twi1.skyframe)*(twi1.skyframe!=0),
                                                      twi1.image/twi1.skyframe, 
                                                                         0.0), 
                                     twi1.header, outname)
                outname = op.join(args.twi_df['Output'][ind], 
                                  'normtwi_%s_%s.fits' 
                                  %(args.twi_df['Specid'][ind],
                                    config.Amp_dict[amp][0]))
                make_amplifier_image(np.where(
                                 np.isfinite(twi2.skyframe)*(twi2.skyframe!=0),
                                                      twi2.image/twi2.skyframe, 
                                                                         0.0), 
                                     twi2.header, outname)
                D = recalculate_dist_coeff(D, twi1, twi2)
                outname2 = op.join(args.twi_df['Output'][ind], 
                                  'mastertrace_%s_%s.dist' 
                                  %(args.twi_df['Specid'][ind],
                                    config.Amp_dict[amp][1]))
                D.writeto(outname2)
                twi1.save_fibers()
                twi2.save_fibers()
                if args.debug:
                    print("Finished working on Cal for %s, %s" %(spec, amp))  

def make_library_image(amp_image, header, outname, fits_list, for_bias=True):
    a,b = amp_image.shape
    overscan = []
    trimsec = []
    bias = re.split('[\[ \] \: \,]', header['BIASSEC'])[1:-1]
    biassec = [int(t)-((i+1)%2) for i,t in enumerate(bias)]  
    trim = re.split('[\[ \] \: \,]', header['TRIMSEC'])[1:-1]
    trimsec = [int(t)-((i+1)%2) for i,t in enumerate(trim)]  
    for F in fits_list:
        overscan.append(biweight_location(F[0].data[biassec[2]:biassec[3],
                                                    biassec[0]:biassec[1]]))
        del F[0].data
    A = []                
    for j,hdu in enumerate(fits_list):
        if for_bias:
            A.append(biweight_filter2d(hdu[0].data[:,i], (25,5),(5,1)) 
            - overscan[j])
        else:
            A.append(hdu[0].data - overscan[j])
    amp_image[:,i] = biweight_location(A, axis=(0,))
    good = np.isfinite(amp_image[:,i])
    amp_image[:,i] = np.interp(np.arange(a), np.arange(a)[good], 
                                   amp_image[good,i])
    for hdu in fits_list:
        del hdu[0].data
        hdu.close()

    hdu = fits.PrimaryHDU(np.array(amp_image[trimsec[2]:trimsec[3],
                                             trimsec[0]:trimsec[1]], 
                                   dtype='float32'), 
                          header=header)
    hdu.header.remove('BIASSEC')
    hdu.header.remove('TRIMSEC')
    hdu.header['DATASEC'] = '[%i:%i,%i:%i]' %(1,trimsec[1]-trimsec[0],1,a)
    hdu.writeto(outname, overwrite=True)  
    
def make_masterbias(args):
    for spec in args.specid:
        spec_ind_bia = np.where(args.bia_df['Specid'] == spec)[0]
        for amp in config.Amps:
            if args.debug:
                print("Working on Bias for %s, %s" %(spec, 
                                                     config.Amp_dict[amp][1]))
                t1 = time.time()
            amp_ind_bia = np.where(args.bia_df['Amp'] == amp)[0]
            bia_sel = np.intersect1d(spec_ind_bia, amp_ind_bia)
            fits_list1 = []
            fits_list2 = []
            for ind in bia_sel:
                fits_list1.append(fits.open(args.bia_df['Files'][ind]))
                fits_list2.append(fits.open(args.bia_df['Files'][ind].replace(amp, 
                                                      config.Amp_dict[amp][0])))
            amp_image = np.zeros((fits_list1[0][0].data.shape))
            header = fits_list1[0][0].header
            outname = op.join(args.configdir, 'lib_bias', 
                              args.bias_outfolder, 'masterbias_%s_%s.fits' 
                                  %(args.bia_df['Specid'][ind], amp)) 
            make_library_image(amp_image, header,outname, fits_list1)
            if args.debug:
                t2 = time.time()
                print('Time Taken for Bias %s, %s: %0.3f' %(spec, amp, t2-t1))
            amp_image = np.zeros((fits_list2[0][0].data.shape))
            header = fits_list2[0][0].header
            outname = op.join(args.configdir, 'lib_bias', 
                              args.bias_outfolder, 'masterbias_%s_%s.fits' 
                                  %(args.bia_df['Specid'][ind],
                                    config.Amp_dict[amp][0])) 
            make_library_image(amp_image, header,outname, fits_list2)
            if args.debug:
                t3 = time.time()
                print('Time Taken for Bias %s, %s: %0.3f' %(spec, 
                                                      config.Amp_dict[amp][0], 
                                                      t3-t2))            

def make_masterdark(args):
    for spec in args.specid:
        spec_ind_drk = np.where(args.drk_df['Specid'] == spec)[0]
        for amp in config.Amps:
            if args.debug:
                print("Working on Dark for %s, %s" %(spec, 
                                                      config.Amp_dict[amp][1])) 
                t1=time.time()                    
            amp_ind_drk = np.where(args.drk_df['Amp'] == amp)[0]
            bia_sel = np.intersect1d(spec_ind_drk, amp_ind_drk)
            fits_list1 = []
            fits_list2 = []
            for ind in bia_sel:
                fits_list1.append(fits.open(args.drk_df['Files'][ind]))
                fits_list2.append(fits.open(args.drk_df['Files'][ind].replace(amp, 
                                                      config.Amp_dict[amp][0])))
            amp_image = np.zeros((fits_list1[0][0].data.shape))
            header = fits_list1[0][0].header
            outname = op.join(args.configdir, 'lib_dark', 
                              args.dark_outfolder, 'masterdark_%s_%s.fits' 
                                  %(args.drk_df['Specid'][ind], amp)) 
            make_library_image(amp_image, header,outname, fits_list1, 
                               for_bias=False)
            if args.debug:
                t2 = time.time()
                print('Time Taken for Dark %s, %s: %0.3f' %(spec, amp, t2-t1))
            amp_image = np.zeros((fits_list2[0][0].data.shape))
            header = fits_list2[0][0].header
            outname = op.join(args.configdir, 'lib_dark', 
                              args.dark_outfolder, 'masterdark_%s_%s.fits' 
                                  %(args.drk_df['Specid'][ind],
                                    config.Amp_dict[amp][0])) 
            make_library_image(amp_image, header,outname, fits_list2, 
                               for_bias=False) 
            if args.debug:
                t3 = time.time()
                print('Time Taken for Dark %s, %s: %0.3f' %(spec, 
                                                      config.Amp_dict[amp][0], 
                                                      t3-t2))  
                                                      
def custom(args):
    lowfib = int(112 / 4. - 1.)
    midfib = int(112 / 2. - 1.)
    highfib = int(3.* 112. / 4. - 1.)
    trace_list = {"LL":[],"LU":[],"RU":[],"RL":[]}
    amps = {"LL":"LL","LU":"LL","RU":"RU","RL":"RU"}
    for spec in args.specid:
        spec_ind_twi = np.where(args.twi_df['Specid'] == spec)[0]
        spec_ind_sci = np.where(args.sci_df['Specid'] == spec)[0]
        for ind in spec_ind_twi:
            amp = args.twi_df['Amp'][ind]
            AMP = amps[amp]
            twi = Amplifier(args.twi_df['Files'][ind],
                            args.twi_df['Output'][ind],
                            calpath=args.twi_df['Output'][ind], 
                            debug=True, dark_mult=0.0,
                            darkpath=args.darkdir, biaspath=args.biasdir,
                            virusconfig=args.configdir, 
                            specname=args.specname[AMP],
                            use_pixelflat=(args.pixelflats<1),
                            init_lims=args.wvl_dict[AMP], 
                            check_fibermodel=True, check_wave=True,
                            fsize=args.fsize, 
                            fibmodel_nbins=args.fibmodel_bins,
                            sigma=args.fibmodel_sig,
                            power=args.fibmodel_pow,
                            use_trace_ref=args.use_trace_ref)
               
            twi.load_fibers()
            if len(twi.fibers)==0:
                twi.get_trace()
            else:
                if not hasattr(twi.fibers[0],'trace'):        
                    twi.get_trace()
            blue = int(twi.D /4.)
            green = int(twi.D /2.)
            red = int(3.*twi.D /4.)
            trace_list[amp].append(np.array([twi.fibers[lowfib].trace[blue],
                                            twi.fibers[lowfib].trace[green],
                                            twi.fibers[lowfib].trace[red],
                                            twi.fibers[midfib].trace[blue],
                                            twi.fibers[midfib].trace[green],
                                            twi.fibers[midfib].trace[red],
                                            twi.fibers[highfib].trace[blue],
                                            twi.fibers[highfib].trace[green],
                                            twi.fibers[highfib].trace[red]]))
        for ind in spec_ind_sci:
            amp = args.sci_df['Amp'][ind]
            AMP = amps[amp]
            print(args.sci_df['Files'][ind])
            
            sci = Amplifier(args.sci_df['Files'][ind],
                            args.sci_df['Output'][ind],
                            calpath=args.twi_dir, skypath=args.sky_dir,
                            debug=False, refit=True, 
                            dark_mult=args.dark_mult[AMP],
                            darkpath=args.darkdir, biaspath=args.biasdir,
                            virusconfig=args.configdir, 
                            specname=args.specname[AMP],
                            use_pixelflat=(args.pixelflats<1),
                            use_trace_ref=args.use_trace_ref,
                            calculate_shift=False)
            sci.load_fibers()
            if len(sci.fibers)==0:
                sci.get_trace()
            else:
                if not hasattr(sci.fibers[0],'trace'):        
                    sci.get_trace()
            blue = int(sci.D /4.)
            green = int(sci.D /2.)
            red = int(3.*sci.D /4.)
            trace_list[amp].append(np.array([sci.fibers[lowfib].trace[blue],
                                            sci.fibers[lowfib].trace[green],
                                            sci.fibers[lowfib].trace[red],
                                            sci.fibers[midfib].trace[blue],
                                            sci.fibers[midfib].trace[green],
                                            sci.fibers[midfib].trace[red],
                                            sci.fibers[highfib].trace[blue],
                                            sci.fibers[highfib].trace[green],
                                            sci.fibers[highfib].trace[red]]))
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12,12))
    ax1 = plt.axes([0.1,0.1,0.35,0.35])
    ax2 = plt.axes([0.1,0.55,0.35,0.35])
    ax3 = plt.axes([0.55,0.1,0.35,0.35])
    ax4 = plt.axes([0.55,0.55,0.35,0.35])
    amps = ["LL","LU","RU","RL"]
    ax = [ax1,ax2,ax3,ax4]
    for i,amp in enumerate(amps):
        TR = np.array(trace_list[amp])
        avg = biweight_location(TR,axis=(0,))
        print(TR-avg)
        ax[i].plot(TR-avg)
    fn = op.join(args.output, args.scidir_date[0], args.instr, 
                        'trace_%s.png' %args.specid[0])
    plt.savefig(fn,dpi=150)
             
def main():
    args = parse_args()
    if args.debug:
        t1 = time.time()
    if args.make_masterbias:
        make_masterbias(args)
    if args.make_masterdark:
        make_masterdark(args)
    if args.custom:
        args.reduce_twi = False
        args.reduce_sci = False
        custom(args)
    if args.reduce_twi:
        reduce_twighlight(args)
    if args.reduce_sci:
        reduce_science(args)                                        
    if args.debug:
        t2=time.time()
        print("Total Time taken: %0.2f s" %(t2-t1))
    

if __name__ == '__main__':
    main()    
