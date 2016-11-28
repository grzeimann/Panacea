# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:17:03 2016

@author: gregz
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from amplifier import Amplifier
from utils import biweight_filter
from fiber_utils import calculate_wavelength_chi2

interactive=True
nbins = 21
smooth_length=21
niter = 1
wavebuff = 100
plotbuff = 70
fiber = 55
#solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/medium_sun.spec')
#gauss = Gaussian1DKernel(13)
#conv = convolve(solar_spec[:,1],gauss)
#solar_spec[:,1] = biweight_filter(conv, int(51/(.1))) / conv
solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/virus_temp.txt')
#solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/lrs2_orange_temp.txt')


orange = Amplifier('/Users/gregz/cure/lrs2_raw/20160806/lrs2/lrs20000001/exp02/lrs2/20160806T020439.7_056RL_twi.fits', 
          '/Users/gregz/cure/lrs2_reductions/twi',
          calpath='/Users/gregz/cure/lrs2_reductions/twi', 
          debug=True, refit=True)
red = Amplifier('/Users/gregz/cure/lrs2_raw/20160806/lrs2/lrs20000001/exp02/lrs2/20160806T020439.7_066LL_twi.fits', 
          '/Users/gregz/cure/lrs2_reductions/twi',
          calpath='/Users/gregz/cure/lrs2_reductions/twi', 
          debug=True, refit=True)
farred = Amplifier('/Users/gregz/cure/lrs2_raw/20160806/lrs2/lrs20000001/exp02/lrs2/20160806T020439.7_066RU_twi.fits', 
          '/Users/gregz/cure/lrs2_reductions/twi',
          calpath='/Users/gregz/cure/lrs2_reductions/twi', 
          debug=True, refit=True)
virus = Amplifier('/Users/gregz/cure/virus_raw/20160512/virus/virus0000001/exp01/virus/20160512T020351.1_084LL_sci.fits',
                  '/Users/gregz/cure/virus_reductions/twi',
                  calpath='/Users/gregz/cure/virus_reductions/twi', 
                  debug=True, refit=True)

instruments = [virus]#, orange, red, farred]
wavelims = [[3480,5500]]#,[4550,7000],[6425,8440],[8230,10500]]
    
for i, instr in enumerate(instruments):
    instr.load_fibers()
    if not instr.fibers:
        instr.fiberextract(use_default_profile=True)
        instr.save_fibers()
    else:
        if instr.fibers[0].spectrum is None:
            instr.fiberextract(use_default_profile=True)
            instr.save_fibers()
    masterwave = []
    masterspec = []
    for fib, fiber in enumerate(instr.fibers):
        if fib==0:
            fiber.wavelength, fiber.wave_polyvals = calculate_wavelength_chi2(np.arange(instr.D), fiber.spectrum, solar_spec, 
                                             init_lims=wavelims[i], 
                                             interactive=interactive, nbins=nbins,smooth_length=smooth_length,
                                             fixscale=False)
        else:
            fiber.wavelength, fiber.wave_polyvals = calculate_wavelength_chi2(np.arange(instr.D), fiber.spectrum, solar_spec, 
                                             init_lims=wavelims[i], 
                                             interactive=False,smooth_length=smooth_length,
                                             init_sol=instr.fibers[fib-1].wave_polyvals)
        y = (biweight_filter(fiber.spectrum, smooth_length) / fiber.spectrum)
        masterwave.append(fiber.wavelength)
        masterspec.append(y)
        plt.scatter(fiber.wavelength, y, color='b', alpha=0.2, edgecolor='none')
    instr.save_fibers()
    masterwave = np.hstack(masterwave)
    masterspec = np.hstack(masterspec)
    ind = np.argsort(masterwave)
    masterwave[:] = masterwave[ind]
    masterspec[:] = masterspec[ind]
    smoothed = biweight_filter(masterspec, 51)
    plt.plot(solar_spec[:,0],solar_spec[:,1],'r-')
    plt.plot(masterwave,smoothed,'k-')
    plt.axis([5300,5500,0.75,1.5])

for fib,fiber in enumerate(virus2.fibers):
    plt.step(fiber.wavelength,fiber.spectrum,color='b',alpha=0.1,where='mid')
plt.axis([3500,4000,0,800])