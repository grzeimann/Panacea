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
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator
import seaborn as sns

# Plot style
sns.set_context('talk')
sns.set_style('ticks')


def get_illum_through(date, fn):
    for gc in ['gc1', 'gc2']:
        t = tarfile.open('/work/03946/hetdex/maverick/%s/gc1/gc1.tar' % date)
        names = t.getnames()
        fns = fn.replace('-', '').replace(':', '')[:-5] + '_gc1_sci.fits'
        N = names + [fns]
        v = np.sort(N)
        ind = np.where(v == fns)[0][0]
        try:
            f = fits.open(t.extractfile(v[ind+1]))
            illum = f[0].header['PUPILLUM']*1.03
            through = f[0].header['TRANSPAR']
            active = f[0].header['GUIDLOOP'] == 'ACTIVE'
        except:
            illum = 0.0
            through = 0.0
            active = False
        t.close()
        if active:
            continue
    
    return illum, through, active
    
plt.figure(figsize=(20, 5))
n = []
fn = sorted(glob.glob('/work/03946/hetdex/maverick/LRS2/STANDARDS/spectrum_20*orange*'))
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
        try:    
            norm = f[0].data['MILLUM'] / 51.4e4 * f[0].data['THROUGHP']
        except:
            print('Could not get header info from reduction for: %s' % f)
            continue
        n.append(g[0].header['OBJECT'][:-6])
        if ('%s' % name) in g[0].header['OBJECT']:
            dt = f.split('_')[1]
            D = datetime.date(int(dt[:4]), int(dt[4:6]), int(dt[6:8]))
            illum, through, active = get_illum_through(dt, g[0].header['DATE'])
            if active:    
                print("Illumination/Throughput for %s is %0.2f" % (f, illum, through))
            else:
                print("No active guider info for: %s" % f)
                illum = 1.0
                through = 1.0
                continue
            dT.append(D)
            d = np.interp(g[0].data[0], wave, flam)
            s.append(biweight(g[0].data[1] * norm / d) / illum / through)
    plt.plot_date(dT, np.array(s), alpha=0.6, ms=10, marker='*')
plt.ylim([0, 1.2])
plt.xlim([datetime.date(2018, 10, 1), datetime.date(2020, 9, 1)])
plt.gcf().autofmt_xdate()

# Plot formatters
myFmt = mdates.DateFormatter('%m/%d/%y')
weeks = mdates.WeekdayLocator()  # every week
months = mdates.MonthLocator()  # every month
days = mdates.DayLocator()  # every day
mL = MultipleLocator(0.1)
ML = MultipleLocator(0.5)

plt.ylabel('Transparency')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(months)
plt.gca().xaxis.set_minor_locator(weeks)
plt.gca().yaxis.set_minor_locator(mL)
plt.gca().yaxis.set_major_locator(ML)
plt.gca().tick_params(axis='x', which='minor', direction='in', bottom=True)
plt.gca().tick_params(axis='x', which='major', direction='in', bottom=True)
plt.gca().tick_params(axis='y', which='minor', direction='in', left=True)
plt.gca().tick_params(axis='y', which='major', direction='in', left=True)
plt.gca().tick_params(axis='both', which='major', length=15, width=3)
plt.gca().tick_params(axis='both', which='minor', length=6, width=1)
plt.savefig('date_throughput.png', dpi=300)