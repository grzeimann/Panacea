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
from fiber_utils import find_maxima

nbins = 20
fiber = 55
solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/medium_sun.spec')


orange = Amplifier('/Users/gregz/cure/lrs2_raw/20160806/lrs2/lrs20000001/exp02/lrs2/20160806T020439.7_056RU_twi.fits', 
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
virus = Amplifier('/Users/gregz/cure/virus_raw/20160512/virus/virus0000001/exp01/virus/20160512T020351.1_075LL_sci.fits',
                  '/Users/gregz/cure/virus_reductions/twi',
                  calpath='/Users/gregz/cure/virus_reductions/twi', 
                  debug=True, refit=True)
       


instruments = [virus]#, orange, red, farred]
wavelims = [[3500,5500]]#,[4550,7000],[6450,8470],[8230,10500]]

def str2bool(v):
  return v.lower() in ("y", "yes", "true", "t", "1")  

for i, instr in enumerate(instruments):
    if not instr.fibers:
        instr.fiberextract()
        instr.save_fibers()
    else:
        if instr.fibers[0].spectrum is None:
            instr.fiberextract()
            instr.save_fibers()
    scale = (wavelims[i][1] - wavelims[i][0])/instr.D
    
    y = (biweight_filter(instr.fibers[fiber].spectrum) 
         / instr.fibers[fiber].spectrum)
    x = np.arange(instr.D)
    peaks, heights = find_maxima(x,y)
    bins = np.linspace(wavelims[i][0], wavelims[i][1], nbins)
    content = False
    for j in xrange(len(bins)):
        while content is False:
            xl = np.searchsorted(solar_spec[:,0],wvi-window_size)
            xu = np.searchsorted(solar_spec[:,0],wvi+window_size)
            ax = plt.axes([0.1,0.1,0.8,0.8])
            ax.step(solar_spec[xl:xu,0],solar_spec[xl:xu,1],where='mid')
            ax.step(wv,yp,'r-',where='mid')
            ax.scatter(wvi, solar_peaks[ind,1])
            #psel = ((solar_peaks[:,0]>wvi-window_size)*
            #        (solar_peaks[:,0]<wvi+window_size))
            #ax.scatter(solar_peaks[psel,0],solar_peaks[psel,1])
            ax.scatter(p_loc_wv[selp],p_height[selp],c='r')
            for s in selp:
                ax.text(p_loc_wv[s],p_height[s]+.03,"%0.2f" %p_loc[s])
            ax.set_xlim([wvi-window_size,wvi+window_size])
            mn = solar_spec[xl:xu,1].min()
            mx = np.max([solar_spec[xl:xu,1].max(),np.max(p_height[selp])])
            rn = mx - mn
            ax.set_ylim([-.1*rn + mn, 1.1*rn+mn])
            plt.show()
            answer = raw_input("Are you happy with the fit?")
            content = str2bool(answer)
            if content:
                continue
            answer = raw_input("What is the xpos for wavelength {:0.0f}?".format(wvi))
            
            plt.close()
            try:
                xv = float(answer)
                matches.append([xv, wvi])
            except ValueError:
                if answer=='off':
                    interactive=False
                elif answer=='q':
                    sys.exit(1)
                else:
                    continue 
gauss = Gaussian1DKernel(13)
conv = convolve(solar_spec[:,1],gauss)
    
# STEP THROUGH THE SOLAR compared with VIRUS, UV, ORANGE, RED, FAR RED for choice of lines
# So I need all 5 to be reduced.

