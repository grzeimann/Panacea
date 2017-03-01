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

import numpy as np
import os.path as op
from astropy.io import fits

from args import parse_args
from amplifier import Amplifier
from utils import biweight_location
import config
import glob

def write_to_fits(hdu, outname):
    try:
        hdu.writeto(outname, overwrite=True)
    except TypeError:
        hdu.writeto(outname, clobber=True)


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
        write_to_fits(hdu, outname)
        hdu = fits.PrimaryHDU(np.array(zcol, dtype='float32'))
        write_to_fits(hdu, outname2)

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
                    sel = w > 1e-4
                    zgrid[k,j,i] = np.sum(data[sel,k]*w[sel])/np.sum(w[sel])
        hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
        zcol = biweight_location(zgrid[:,:,:],axis=(0,))
        hdu.header['CDELT3'] = F1[0].header['CDELT1']
        hdu.header['CRVAL3'] = F1[0].header['CRVAL1']
        hdu.header['CRPIX3'] = F1[0].header['CRPIX1']
        write_to_fits(hdu, outname)
        hdu = fits.PrimaryHDU(np.array(zcol, dtype='float32'))
        write_to_fits(hdu, outname2)
        
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
    write_to_fits(hdu, outname)


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
                              usecols=[0,1,2,4], skiprows=args.ifucen_fn[amp][1])
                if args.debug:
                    print("Working on Sci/Twi for %s, %s" %(spec, amp)) 
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
                                 cont_smooth=args.cont_smooth,
                                 make_residual=False, do_cont_sub=False,
                                 make_skyframe=False)
                sci1.load_all_cal()
                wavelim=[4500,4600]
                xlim = np.interp([wavelim[0],wavelim[1]],
                                 np.linspace(args.wvl_dict[amp][0],
                                             args.wvl_dict[amp][1], sci1.D),
                                 np.arange(sci1.D))
                cols=np.arange(int(xlim[0])-10,int(xlim[1])+10)  
                sci1.fiberextract(cols=cols)
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
                                 cont_smooth=args.cont_smooth,
                                 make_residual=False, do_cont_sub=False,
                                 make_skyframe=False)
                sci2.load_all_cal()
                sci2.fiberextract(cols=cols)                    
                sci2.sky_subtraction()
                Fe, FeS = recreate_fiberextract(sci1, sci2, 
                                                wavelim=wavelim, 
                                                disp=args.disp[amp])
                FE = [Fe, FeS]
                FEN = ['Fe', 'FeS']
                for f,n in zip(FE, FEN):
                    outname = op.join(args.sci_df['Output'][ind],
                                      '%s%s_%s_sci_%s.fits' %(n,
                              op.basename(args.sci_df['Files'][ind]).split('_')[0],
                                                       args.sci_df['Ifuslot'][ind], 
                                                          config.Amp_dict[amp][1]))
                    make_fiber_image(f, sci1.header, outname, args, amp)
                    
                    make_cube_file(args, outname, ifucen, args.cube_scale, 
                                   config.Amp_dict[amp][1])
                if args.save_sci_fibers:
                    sci1.save_fibers()
                    sci2.save_fibers()
                if args.save_sci_amplifier:
                    sci1.save()
                    sci2.save()
                if args.debug:
                    print("Finished working on Sci/Twi for %s, %s" %(spec, amp))        

def main():
    args = parse_args()
    if args.debug:
        t1 = time.time()
    if args.reduce_sci:
        reduce_science(args)                                        
    if args.debug:
        t2 = time.time()
        print("Total Time taken: %0.2f s" %(t2-t1))
    

if __name__ == '__main__':
    main()    