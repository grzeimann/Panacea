#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 14:48:58 2022

@author: gregz
"""

import numpy as np
from astropy.io import fits
import glob
import os.path as op
from astropy.coordinates import get_moon, get_sun
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import EarthLocation

def moon_phase_angle(time, location=None, ephemeris=None):
    """
    Calculate lunar orbital phase in radians.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        `~astropy.coordinates.solar_system_ephemeris` (which is
        set to 'builtin' by default).

    Returns
    -------
    i : `~astropy.units.Quantity`
        Phase angle of the moon [radians]
    """

    sun = get_sun(time)
    moon = get_moon(time, location=location, ephemeris=ephemeris)
    elongation = sun.separation(moon)
    return np.arctan2(sun.distance*np.sin(elongation),
                      moon.distance - sun.distance*np.cos(elongation)), moon



def moon_illumination(time, location=None, ephemeris=None):
    """
    Calculate fraction of the moon illuminated.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        `~astropy.coordinates.solar_system_ephemeris` (which is
        set to 'builtin' by default).

    Returns
    -------
    k : float
        Fraction of moon illuminated
    """
    i, moon = moon_phase_angle(time, location=location, ephemeris=ephemeris)
    k = (1 + np.cos(i))/2.0
    return k.value, moon, i


filenames = np.array(sorted(glob.glob('/work/03946/hetdex/maverick/LRS2/UT22-1-002/multi*uv.fits')))
dates = [op.basename(f).split('_')[1] for f in filenames]
dates = np.array(dates)

wave = np.arange(3650, 10500, 0.7)
skies = []
dateobs = []
ras = []
decs = []
loc = EarthLocation.of_site('McDonald Observatory')
for filename in filenames:
    f = fits.open(filename)
    name = f[0].header['OBJECT']
    millum = f[0].header['MILLUM']
    throughp = f[0].header['THROUGHP']
    
    if millum == 51e4:
        continue
    if (throughp < 0.1) + (throughp == 1.):
        continue
    slot = name.split('_')[-2]
    f = fits.open(filename.replace('uv', 'red'))
    g = fits.open(filename.replace('uv', 'farred'))
    skyred = f[1].data[125]
    skyfarred = g[1].data[125]
    wavered = f[6].data[0]
    wavefarred = g[6].data[0]
    skyre = np.interp(wave, wavered, skyred, 
                      left=np.nan, right=np.nan)
    skyfr = np.interp(wave, wavefarred, skyfarred, 
                      left=np.nan, right=np.nan)
    f = fits.open(filename.replace('uv', 'uv'))
    g = fits.open(filename.replace('uv', 'orange'))
    skyuv = f[1].data[125]
    skyorange = g[1].data[125]
    waveuv = f[6].data[0]
    waveorange = g[6].data[0]
    skyuv = np.interp(wave, waveuv, skyuv, 
                      left=np.nan, right=np.nan)
    skyor = np.interp(wave, waveorange, skyorange, 
                      left=np.nan, right=np.nan)
    wsel = (wave < 6880) * (wave>6520)
    norm = np.nanmedian(skyre[wsel] / skyor[wsel])
    print('Norm for %s %s: %0.2f' % (filename, name, norm))
    skyuv *= norm
    skyor *= norm
    sky = np.nanmean([skyuv, skyor, skyre, skyfr], axis=0)
    if np.abs(norm - 1.) < 0.5:
        skies.append(sky)
        ras.append(f[0].header['QRA'])
        decs.append(f[0].header['QDEC'])
        dateobs.append(f[0].header['DATE'])
skies = np.array(skies)
fits.HDUList([fits.PrimaryHDU(skies), fits.BinTableHDU(Table([dateobs, ras, decs], names=['Date', 'RA', 'Dec']))]).writeto('skyfile3.fits', overwrite=True)
