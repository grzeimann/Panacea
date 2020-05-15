#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:53:16 2020

@author: gregz
"""
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import sys
from astropy.convolution import convolve, Gaussian1DKernel
from matplotlib.ticker import MultipleLocator
from astropy.table import Table
from scipy.signal import savgol_filter


filename = sys.argv[1]
filename2 = sys.argv[2]
# Plot Style
sns.set_context('talk')
sns.set_style('ticks')
plt.figure(figsize=(10, 6))
nexps = 1
waves_dict = {'uv': np.array([3670., 4050., 4200., 4420., 4500., 4600.]),
              'orange': np.array([4750., 5050., 5600., 6000., 6650.]),
              'red': np.array([6700., 7100., 7450., 7700., 8000.]),
              'farred': np.array([8450, 8600., 8800., 9900., 10100.])}
def get_cor(calbase, channel):
    waves = waves_dict[channel]
    k = fits.open('%s_exp01_%s.fits' % (calbase, channel))
    if 'red' in channel:
        objname = k[0].header['OBJECT'].split('_066_')[0].lower()
    else:
        objname = k[0].header['OBJECT'].split('_056_')[0].lower()
    m = np.loadtxt('/work/03946/hetdex/maverick/virus_config/standards/m%s.dat.txt' % objname)
    flux = 10**(-0.4 * (m[:, 1]-23.9))*1e-29 * 2.99792e18 / m[:, 0]**2
    if objname == 'gd248':
        sel = m[:, 0] > 6250.
        s = savgol_filter(flux[sel], 425, 2)
        t = np.hstack([flux[~sel], s])
        plt.plot(m[:, 0], t)
        sel = m[:, 0] > 8500.
        p = np.polyfit(m[sel, 0], 1./t[sel], 3)
        d = m[-1, 0] - m[-2, 0]
        x = np.arange(m[-1, 0]+d, 10500+d, d)
        M = np.hstack([m[:, 0], x])
        Y = np.hstack([t, 1./np.polyval(p, x)])
        model = np.interp(k[0].data[0], M, Y)
    else:
        model = np.interp(k[0].data[0], m[:, 0], flux)
    p = np.polyfit(waves, np.interp(waves, k[0].data[0], k[0].data[1] / model), 2)
    init = np.polyval(p, k[0].data[0])
    cor = k[0].data[1] / (model * init)
    cor[np.abs(k[0].data[0]-6563.)<30.] = 1.
    cor[np.abs(k[0].data[0]-4861.)<30.] = 1.
    cor[k[0].data[0] < 4600.] = 1.
    return init, cor

def_wave = np.arange(3650., 10500., 0.7)
spec = {'uv': [], 'orange': [], 'red': [], 'farred': []}
wave = {'uv': None, 'orange': None, 'red': None, 'farred': None}
Spec = []
Cor = []
Err = []
Sky = []
allspec = []
allerr = []
allsky = []
c = []
for base, calbase, nexp, channels in zip([filename],
                               [filename2],
                               [1],
                               [['uv', 'orange']]):
    
    for channel in channels:
        cor, CO = get_cor(calbase, channel)
        for i in np.arange(1, nexp+1):
            f = fits.open('%s_exp%02d_%s.fits' % (base, i, channel))
            t = np.interp(def_wave, f[0].data[0], f[0].data[1] / CO / cor, left=0., right=0.)
            s = np.interp(def_wave, f[0].data[0], f[0].data[2] / CO / cor, left=0., right=0.)
            e = np.interp(def_wave, f[0].data[0], f[0].data[-1] / CO / cor, left=0., right=0.)
            spec[channel].append(t)
            wave[channel] = def_wave
            allspec.append(t)
            allerr.append(e)
            allsky.append(s)
        c.append(np.interp(def_wave, f[0].data[0], CO*cor, left=0., right=0.))
allspec = np.array(allspec)
allspec[allspec==0.] = np.nan
#
#lims = [[6450., 6950., 6450., 7000.], [8275., 8400., 8100., 8550.]]
#for j, a in enumerate(allspec):
#    for lim in lims:
#        sel = np.isfinite(a) * (def_wave>lim[0]) * (def_wave<lim[1])
#        if sel.sum():
#            left = np.nanmedian(allspec[:, np.abs(def_wave-lim[2])<100.])
#            right = np.nanmedian(allspec[:, np.abs(def_wave-lim[3])<100.])
#            print(left, right, lim, j)
#            m = (right - left) / (lim[3] - lim[2])
#            d = np.polyval(np.polyfit(def_wave[sel], a[sel], 2), def_wave[sel])
#            y = m * (def_wave[sel] - lim[2]) + left
#            allspec[j][sel] = a[sel] * y / d
allerr = np.array(allerr)
allerr[allerr==0.] = np.nan
allsky = np.array(allsky)
allsky[allerr==0.] = np.nan
c = np.array(c)
c[c==0.] = np.nan
Spec.append(np.nanmean(allspec, axis=0))
Err.append(np.nanmean(allerr, axis=0) / np.sqrt(np.isfinite(allerr).sum(axis=0)))
Sky.append(np.nanmean(allsky, axis=0))
Cor.append(np.nanmean(c, axis=0))
plt.plot(def_wave, convolve(np.nanmean(allspec, axis=0), Gaussian1DKernel(2.8)), lw=1.0, alpha=0.4, label=base, zorder=2)
Table([def_wave, Spec[-1], Err[-1], Sky[-1], Cor[-1]], names=['wavelength', 'f_lam', 'e_lam', 'sky_lam', 'tel_cor']).write(base+'_coadd.txt', overwrite=True, format='ascii.fixed_width_two_line')
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.gca().tick_params(axis='y', which='both', left=True, right=True)
plt.gca().tick_params(axis='x', which='both', bottom=True, top=True)
plt.gca().tick_params(axis='both', which='major', length=15, width=3)
plt.gca().tick_params(axis='both', which='minor', length=6, width=2)
ML = MultipleLocator(1000)
ml = MultipleLocator(200)
ran = (np.nanpercentile(Spec[-1], 98) - np.nanpercentile(Spec[-1], 2))
low = np.nanpercentile(Spec[-1], 98)
plt.gca().xaxis.set_major_locator(ML)
plt.gca().xaxis.set_minor_locator(ml)
plt.xlabel('Wavelength')
plt.ylabel('F$_{\lambda}$')
plt.legend()
plt.xlim([3650, 7000])
plt.ylim([low-ran*0.2, low+ran*1.2])
plt.savefig('%s.png' % filename, dpi=300)