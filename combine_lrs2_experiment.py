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
from scipy.signal import savgol_filter, medfilt
from math_utils import biweight


filenames = [v.replace(' ', '') for v in sys.argv[1].split(',')]
filenames2 = [v.replace(' ', '') for v in sys.argv[2].split(',')]
sides = [v.replace(' ', '') for v in sys.argv[3].split(',')]
outfile = sys.argv[4]
sidedict = {'blue': ['uv', 'orange'], 'red': ['red', 'farred']}

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
    cor[np.abs(k[0].data[0]-6563.)<50.] = np.nan
    cor[np.abs(k[0].data[0]-4861.)<50.] = np.nan
    cor[np.isnan(cor)] = np.interp(k[0].data[0][np.isnan(cor)], k[0].data[0][np.isfinite(cor)],
                                   cor[np.isfinite(cor)], left=1.0, right=1.0)
    cor[k[0].data[0] < 4550.] = 1.
    return init, cor

def connect_channels(spec1, spec2, def_wave, w1, w2, w3, lw, hw):
    niter = 3
    for j in np.arange(niter):
        n3 = np.nanmedian(spec1[np.abs(def_wave-w1)<10.])
        n4 = np.nanmedian(spec2[np.abs(def_wave-w2)<10.])
        n5 = np.nanmedian(spec2[np.abs(def_wave-w3)<10.])
        sel = (def_wave > lw) * (def_wave < hw)
        sel1 = sel * np.isfinite(spec1)
        sel2 = sel * np.isfinite(spec2)
        p0 = np.polyfit([w1, w2, w3], [n3, n4, n5], 2)
        p1 = np.polyfit(def_wave[sel1], spec1[sel1], 2)
        p2 = np.polyfit(def_wave[sel2], spec2[sel2], 2)
        norm = np.polyval(p0, def_wave[sel])
        norm1 = np.polyval(p1, def_wave[sel])
        norm2 = np.polyval(p2, def_wave[sel])
        spec1[sel] = spec1[sel] / norm1 * norm
        spec2[sel] = spec2[sel] / norm2 * norm
        nl = np.nanmedian(spec1[np.abs(def_wave-(lw-3.))<3.])
        nh = np.nanmedian(spec1[np.abs(def_wave-(lw+3.))<3.])
        mult = nh / nl / 1.01
        spec1[def_wave<=lw] = spec1[def_wave<=lw] * mult
        nl = np.nanmedian(spec2[np.abs(def_wave-(hw-3.))<3.])
        nh = np.nanmedian(spec2[np.abs(def_wave-(hw-3.))<3.])
        mult = nl / nh / 1.01
        spec2[def_wave>=hw] = spec2[def_wave>=hw] * mult
    return spec1, spec2

def_wave = np.arange(3650., 10500., 0.7)
wave = {'uv': None, 'orange': None, 'red': None, 'farred': None}
normdict = {'blue': [4260, 4800, 5100, 4580, 4690], 
            'red': [8000, 8600, 8700, 8150, 8560]}
Spec = []
Cor = []
Err = []
Sky = []
allspec = []
redspec = []
bluespec = []
blueerr = []
rederr = []
bluesky = []
redsky = []
c = []
for base, calbase, side in zip(filenames, filenames2, sides):
    channels = sidedict[side]
    nexp = len(glob.glob('%s_exp*_%s.fits' % (base, channels[0])))
    print('Found %i exposures' % nexp)
    for channel in channels:
        cor, CO = get_cor(calbase, channel)
        for i in np.arange(1, nexp+1):
            f = fits.open('%s_exp%02d_%s.fits' % (base, i, channel))
            t = np.interp(def_wave, f[0].data[0], f[0].data[1] / CO / cor, left=0., right=0.)
            s = np.interp(def_wave, f[0].data[0], f[0].data[2] / CO / cor, left=0., right=0.)
            e = np.interp(def_wave, f[0].data[0], f[0].data[-1] / CO / cor, left=0., right=0.)
            wave[channel] = def_wave
            if side == 'blue':
                bluespec.append(t)
                blueerr.append(e)
                bluesky.append(s)
            else:
                redspec.append(t)
                rederr.append(e)
                redsky.append(s)
        c.append(np.interp(def_wave, f[0].data[0], CO*cor, left=0., right=0.))
    N = nexp * len(channels)
    w1, w2, w3, lw, hw = normdict[side]
#    for i in np.arange(nexp):
#        ind1 = -N + i
#        ind2 = -N + nexp + i
#        allspec[ind1], allspec[ind2] = connect_channels(allspec[ind1],
#                                                        allspec[ind2], 
#                                                        def_wave, w1, w2, w3,
#                                                        lw, hw)
    
bluespec = np.array(bluespec)
bluespec[bluespec==0.] = np.nan
redspec = np.array(redspec)
redspec[redspec==0.] = np.nan
blueerr = np.array(blueerr)
blueerr[blueerr==0.] = np.nan
rederr = np.array(rederr)
rederr[rederr==0.] = np.nan
bluesky = np.array(bluesky)
bluesky[bluesky==0.] = np.nan
redsky = np.array(redsky)
redsky[redsky==0.] = np.nan
c = np.array(c)
c[c==0.] = np.nan
Blue = np.nanmean(bluespec, axis=0)
norm = np.nanmedian(bluespec / Blue[np.newaxis, :], axis=1)
bluespec = bluespec / norm[:, np.newaxis]
blueerr = blueerr / norm[:, np.newaxis]
weights = 1. / blueerr**2
weights = weights / np.nansum(weights, axis=0)[np.newaxis, :]
Blue = np.nansum(bluespec*weights, axis=0)
BlueErr = np.sqrt(np.nansum(blueerr**2*weights, axis=0))
BlueSky = np.nansum(bluesky*weights, axis=0)
Blue[np.abs(def_wave-3735.7)<0.5] = np.nan
Blue[Blue==0.] = np.nan
BlueSky[BlueSky==0.] = np.nan
BlueErr[BlueErr==0.] = np.nan

Red = np.nanmean(redspec, axis=0)
norm = np.nanmedian(redspec / Red[np.newaxis, :], axis=1)
redspec = redspec / norm[:, np.newaxis]
rederr = rederr / norm[:, np.newaxis]
weights = 1. / rederr**2
weights = weights / np.nansum(weights, axis=0)[np.newaxis, :]
Red = np.nansum(redspec*weights, axis=0)
RedErr = np.sqrt(np.nansum(rederr**2*weights, axis=0))
RedSky = np.nansum(redsky*weights, axis=0)
Red[Red==0.] = np.nan
RedSky[RedSky==0.] = np.nan
RedErr[RedErr==0.] = np.nan
sel = np.isfinite(Red) * np.isfinite(Blue)
Norm = biweight(Blue[sel] / Red[sel])
print(Norm)
Red *= Norm
RedErr *= Norm
RedSky *= Norm

Spec = np.nanmean([Blue, Red], axis=0)
Err = np.sqrt(np.nansum([BlueErr**2, RedErr**2], axis=0))
Sky = np.nanmean([BlueSky, RedSky], axis=0)
Cor = np.nanmean(c, axis=0)

for spec in [bluespec, redspec]:
    for s in spec:
        plt.plot(def_wave, s, lw=1.0, alpha=0.4, zorder=1)
plt.plot(def_wave, Spec, 'k-', lw=1.0, alpha=0.4, zorder=2)
Table([def_wave, Spec, Err, Sky, Cor], names=['wavelength', 'f_lam', 'e_lam', 'sky_lam', 'tel_cor']).write(outfile+'_coadd.txt', overwrite=True, format='ascii.fixed_width_two_line')
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
locu = np.where(np.isfinite(Spec))[0][-1]
locl = np.where(np.isfinite(Spec))[0][0]
plt.xlim([def_wave[locl], def_wave[locu]])
plt.ylim([low-ran*0.2, low+ran*1.2])
plt.savefig('%s.png' % outfile, dpi=300)