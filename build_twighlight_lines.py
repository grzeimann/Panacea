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
from fiber_utils import calculate_wavelength_chi2, calculate_wavelength

interactive=True
nbins = 25
smooth_length=21
niter = 1
wavebuff = 100
plotbuff = 70
fiber = 55
#solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/medium_sun.spec')
#gauss = Gaussian1DKernel(15)
#conv = convolve(solar_spec[:,1],gauss)
#solar_spec[:,1] = (biweight_filter(conv, int(51/(.1))) / conv - 1)/1.0 + 1
solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/lrs2_uv_temp.txt')
#solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/lrs2_orange_temp.txt')

uv = Amplifier('/Users/gregz/cure/lrs2_raw/20161213/lrs2/lrs20000001/exp01/lrs2/20161213T001016.5_056LU_twi.fits', 
          '/Users/gregz/cure/lrs2_reductions/20161213/lrs2/lrs20000001/exp01/lrs2',
          calpath='/Users/gregz/cure/lrs2_reductions/20161213/lrs2/lrs20000001/exp01/lrs2', 
          debug=True, refit=True, specname="lrs2_uv")
orange = Amplifier('/Users/gregz/cure/lrs2_raw/20161213/lrs2/lrs20000001/exp01/lrs2/20161213T001016.5_056RL_twi.fits', 
          '/Users/gregz/cure/lrs2_reductions/20161213/lrs2/lrs20000001/exp01/lrs2/',
          calpath='/Users/gregz/cure/lrs2_reductions/20161213/lrs2/lrs20000001/exp01/lrs2/', 
          debug=True, refit=True, specname="lrs2_orange")
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

instruments = [uv]#, orange, red, farred]
wavelims = [[3633,4655]]#,[4550,7000],[6425,8440],[8230,10500]]
    
for i, instr in enumerate(instruments):
    #instr.load_fibers()
    if not instr.fibers:
        instr.fiberextract(use_default_profile=False, check_fibermodel=True, 
                           fsize=6., bins=11)
        instr.save_fibers()
    else:
        if instr.fibers[0].spectrum is None:
            instr.fiberextract(use_default_profile=False, check_fibermodel=True)
            instr.save_fibers()
    masterwave = []
    masterspec = []
    plt.figure(figsize=(12,10))
    instr.get_wavelength_solution(check_wave=True, interactive=True, 
                                  init_lims=wavelims[i], filt_size_sky=151, 
                                  nbins=31)
    for fib, fiber in enumerate(instr.fibers):
        y = (biweight_filter(fiber.spectrum, smooth_length) / fiber.spectrum)
        masterwave.append(fiber.wavelength)
        masterspec.append(y)
        plt.scatter(fiber.wavelength, y, color='b', alpha=0.2, edgecolor='none')
    instr.get_fiber_to_fiber()
    instr.save_fibers()
    masterwave = np.hstack(masterwave)
    masterspec = np.hstack(masterspec)
    ind = np.argsort(masterwave)
    masterwave[:] = masterwave[ind]
    masterspec[:] = masterspec[ind]
    smoothed = biweight_filter(masterspec, 151)
    plt.plot(solar_spec[:,0],solar_spec[:,1],'r-')
    plt.plot(masterwave,smoothed,'k-')
    plt.axis([wavelims[i][1]-200,wavelims[i][1]+0,0.75,2.0])

