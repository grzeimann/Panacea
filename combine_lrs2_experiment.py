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
import glob
from astropy.convolution import convolve, Gaussian1DKernel
from matplotlib.ticker import MultipleLocator
from astropy.table import Table
from scipy.signal import savgol_filter
from math_utils import biweight


filename = sys.argv[1]
filename2 = sys.argv[2]
# Plot Style
sns.set_context('talk')
sns.set_style('ticks')
plt.figure(figsize=(10, 6))
nexps = 1
waves_dict = {'uv': np.array([3670., 4050., 4200., 4420., 4500., 4560.]),
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
    Z = waves*0.
    for j, wave in enumerate(waves):
        Z[j] = np.nanmedian((k[0].data[1] / model)[np.abs(k[0].data[0] - wave)<5.])
    p = np.polyfit(waves, Z, 2)
    init = np.polyval(p, k[0].data[0])
    cor = k[0].data[1] / (model * init)
    cor[np.abs(k[0].data[0]-6563.)<30.] = 1.
    cor[np.abs(k[0].data[0]-4861.)<30.] = 1.
    cor[k[0].data[0] < 4550.] = 1.
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
for base, calbase, channels in zip([filename],
                               [filename2],
                               [['uv', 'orange']]):
    nexp = len(glob.glob('%s_exp*_%s.fits' % (base, channels[0])))
    print('Found %i exposures' % nexp)
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

for i in np.arange(nexp):
    a1 = allspec[i]
    a2 = allspec[i+nexp]
    n1 = np.nanmedian(a1[np.abs(def_wave-4640.)<5.])
    n2 = np.nanmedian(a2[np.abs(def_wave-4640.)<5.])
    avg = (n1 + n2) / 2.
    n3 = np.nanmedian(a1[np.abs(def_wave-4610.)<5.])
    n4 = np.nanmedian(a2[np.abs(def_wave-4690.)<5.])
    sel = (def_wave > 4610.) * (def_wave < 4690.)
    sel1 = sel * np.isfinite(a1)
    sel2 = sel * np.isfinite(a2)
    p0 = np.polyfit([4610., 4640., 4690.], [n3, avg, n4], 2)
    p1 = np.polyfit(def_wave[sel1], a1[sel1], 2)
    p2 = np.polyfit(def_wave[sel2], a2[sel2], 2)
    norm = np.polyval(p0, def_wave[sel])
    norm1 = np.polyval(p1, def_wave[sel])
    norm2 = np.polyval(p2, def_wave[sel])
    allspec[i][sel] = allspec[i][sel] / norm1 * norm
    allspec[i+nexp][sel] = allspec[i+nexp][sel] / norm2 * norm
allerr = np.array(allerr)
allerr[allerr==0.] = np.nan
allsky = np.array(allsky)
allsky[allerr==0.] = np.nan
c = np.array(c)
c[c==0.] = np.nan
Spec = np.nanmean(allspec, axis=0)
norm = np.nanmedian(allspec / Spec[np.newaxis, :], axis=1)
allspec = allspec / norm[:, np.newaxis]
allerr = allerr / norm[:, np.newaxis]
weights = 1. / allerr**2
weights = weights / np.nansum(weights, axis=0)[np.newaxis, :]
Spec = np.sum(allspec*weights, axis=0)
Err = np.nanmean(allerr, axis=0) / np.sqrt(np.isfinite(allerr).sum(axis=0))
Sky = np.nanmean(allsky, axis=0)
Cor = np.nanmean(c, axis=0)
Spec[np.abs(def_wave-3735.7)<0.5] = np.nan
Spec[np.abs(def_wave-4620.)<70.] = np.nan

for s in allspec:
    plt.plot(def_wave, s, lw=1.0, alpha=0.4, zorder=1)
plt.plot(def_wave, Spec, 'k-', lw=1.0, alpha=0.4, zorder=2)
Table([def_wave, Spec, Err, Sky, Cor], names=['wavelength', 'f_lam', 'e_lam', 'sky_lam', 'tel_cor']).write(base+'_coadd.txt', overwrite=True, format='ascii.fixed_width_two_line')
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.gca().tick_params(axis='y', which='both', left=True, right=True)
plt.gca().tick_params(axis='x', which='both', bottom=True, top=True)
plt.gca().tick_params(axis='both', which='major', length=15, width=3)
plt.gca().tick_params(axis='both', which='minor', length=6, width=2)
ML = MultipleLocator(1000)
ml = MultipleLocator(200)
ran = (np.nanpercentile(Spec, 98) - np.nanpercentile(Spec, 2))
low = np.nanpercentile(Spec, 2)
plt.gca().xaxis.set_major_locator(ML)
plt.gca().xaxis.set_minor_locator(ml)
plt.xlabel('Wavelength')
plt.ylabel('F$_{\lambda}$')
plt.xlim([3650, 7000])
plt.ylim([low-ran*0.2, low+ran*1.2])
plt.savefig('%s.png' % filename, dpi=300)