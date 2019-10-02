#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:16:51 2019

@author: gregz
"""
import matplotlib
matplotlib.use('agg')
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import datetime
from math_utils import biweight
import tarfile


def get_illum(date, fn):
    t = tarfile('/work/03946/hetdex/maverick/%s/gc1/gc1.tar')
    names = t.getnames()
    fns = fn.replace('-', '').replace(':')[:-5] + '_gc1_sci.fits'
    N = names + [fns]
    v = np.array(sorted(N))
    ind = np.where(v == fns)[0]
    try:
        f = fits.open(t.extractfile(v[ind+1]))
        illum = f[0].header['PUPILLUM']
    except:
        illum = 0.0
    t.close()
    return illum
    
plt.figure(figsize=(10, 5))
n = []
fn = sorted(glob.glob('/work/03946/hetdex/maverick/LRS2/STANDARDS/spectrum_201*orange*'))
for name in ['BD+40_4032', 'BD_+17_4708', 'FEIGE_110', 'FEIGE_34',
             'FEIGE_66', 'FEIGE_67', 'G191B2B', 'G193-74',
             'GD140', 'GD248', 'GD50', 'HD_19445', 'HD_84937', 'HILTNER_102',
             'HILTNER_600', 'HZ_21', 'HZ_4', 'PG0823+546', 'WOLF_1346']:
    
    try:
        T = np.loadtxt('/work/03946/hetdex/maverick/virus_config/standards/m%s.dat.txt' % name.lower())
    except:
        print("Can't load %s" % name)
        continue
    flam = 10**(-0.4*(T[:, 1]-23.9))*1e-29 *3e18 / T[:, 0]**2
    wave = T[:, 0]
    s = []
    dT = []
    for f in fn:
        g = fits.open(f)
        n.append(g[0].header['OBJECT'][:-6])
        if ('%s' % name) in g[0].header['OBJECT']:
            dt = f.split('_')[1]
            D = datetime.date(int(dt[:4]), int(dt[4:6]), int(dt[6:8]))
            if D > datetime.date(2019, 4, 1):
                illum = get_illum(dt, g[0].header['DATE'])
                print("Illumination for %s is %0.2f" % (f, illum))
            dT.append(D)
            d = np.interp(g[0].data[0], wave, flam)
            s.append(biweight(g[0].data[1] / d) / illum)
    plt.plot_date(dT, np.array(s), alpha=0.6, ms=5)
plt.ylim([0, 1.2])
plt.xlim([datetime.date(2018, 10, 1), datetime.date(2019, 10, 1)])
plt.savefig('date_throughput.png', dpi=300)