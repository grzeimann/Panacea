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
import os.path as op
from scipy.ndimage.filters import percentile_filter
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=False)
# Plot style
sns.set_context('talk')
sns.set_style('ticks')


def get_illum_through(date, fn):
    for gc in ['gc1']:
        t = tarfile.open('/work/03946/hetdex/maverick/%s/%s/%s.tar' % (date, gc, gc))
        onames = t.getnames()
        names = [op.basename(name) for name in onames]
        fns = fn.replace('-', '').replace(':', '')[:-5] + ('_%s_sci.fits' % gc)
        N = names + [fns]
        v = np.sort(N)
        ind = np.where(v == fns)[0][0]
        ind = np.where([v[ind+1] == op.basename(name) for name in onames])[0][0]
        try:
            f = fits.open(t.extractfile(onames[ind]))
            illum = f[0].header['PUPILLUM']*1.03
        except:
            illum = 0.0
    return illum
    
plt.figure(figsize=(20, 5))
n = []
fn = sorted(glob.glob('/work/03946/hetdex/maverick/LRS2/STANDARDS/spectrum_*orange*.fits'))
names = ['BD+40_4032', 'BD_+17_4708', 'FEIGE_110', 'FEIGE_34',
             'FEIGE_66', 'FEIGE_67', 'G191B2B', 'G193-74',
             'GD140', 'GD248', 'GD50', 'HD_19445', 'HD_84937', 'HILTNER_102',
             'HILTNER_600', 'HZ_21', 'HZ_4', 'PG0823+546', 'WOLF_1346']
names = ['BD+40_4032', 'BD_+17_4708', 'FEIGE_110', 'FEIGE_34',
             'FEIGE_66', 'FEIGE_67', 'G191B2B', 'G193-74',
             'GD140', 'GD248', 'GD50', 'HD_19445', 'HD_84937', 'HILTNER_102',
             'HILTNER_600', 'HZ_21', 'HZ_4', 'PG0823+546', 'WOLF_1346',
             'BD_+26_2606', 'FEIGE34', 'G193_74', 'HZ_44', 'BD+28_4211',
             'BD+25_3941', 'BD+40_4032', 'BD+33_2642']
cmap = matplotlib.cm.get_cmap('magma')
colors = cmap(np.linspace(0, 1, len(names)))
alldT = []
alls = []
allss = []
for name, color in zip(names, colors):
    
    try:
        T = np.loadtxt('/work/03946/hetdex/maverick/virus_config/standards/m%s.dat.txt' % name.lower())
    except:
        print("Can't load %s" % name)
        continue
    flam = 10**(-0.4*(T[:, 1]-23.9))*1e-29 *3e18 / T[:, 0]**2
    wave = T[:, 0]
    s = []
    ss = []
    dT = []
    for f in fn:
        g = fits.open(f)
        n.append(g[0].header['OBJECT'][:-6])
        if ('%s' % name) in g[0].header['OBJECT']:
            dt = f.split('_')[1]
            try:    
                A = g[0].header['MILLUM'] / 51.4e4
                if np.abs(A - 1.) < 0.01:
                    try:
                        illum = get_illum_through(dt, g[0].header['DATE'])
                        if illum <= 0.0:
                            print('Could not get guider info for %s' % f)
                            continue
                    except:
                        print('Could not get guider info for %s' % f)
                        continue
                    norm = 1. / illum
                else:
                    norm = 1.
            except:
                print('Could not get header info from reduction for: %s' % f)
                continue
            thr_flag = True
            try:
                thr = g[0].header['THROUGHP']
                if (thr < 0.1) + (thr > 1.5):
                    thr_flag = False
                    thr = 0.0
                if thr == 1.:
                    continue
                print('Header throughput for %s: %0.2f' % (f, thr))
                norm *= thr
            except:
                continue
            D = datetime.date(int(dt[:4]), int(dt[4:6]), int(dt[6:8]))
#            try:
#                illum, through, active = get_illum_through(dt, g[0].header['DATE'])
#            except:
#                print('Could not get guider info for %s' % f)
#                continue
#            if (through < 0.1) + (through > 1.5):
#                through = 0.0
#            if illum < 0.1:
#                illum = 1.0
#            print("Illumination/Throughput for %s is %0.2f, %0.2f" % (f, illum, through))
            if not thr_flag:
                continue
            dT.append(D)
            alldT.append(D)
            d = np.interp(g[0].data[0], wave, flam)
            s.append(biweight(g[0].data[1][300:800] * norm / d[300:800]))
            alls.append(s[-1])
            ss.append(thr)
            allss.append(thr)
    plt.plot_date(dT, np.array(s), alpha=0.8, ms=10, marker='*', color='firebrick')
    plt.plot_date(dT, np.array(ss), alpha=0.8, ms=3, marker='s', color='dimgray')
inds = np.argsort(alldT)
S = np.array(alls)[inds]
SS = np.array(allss)[inds]
plt.plot_date(np.array(alldT)[inds], percentile_filter(S, 75, size=50), 'r-', color='tomato', lw=3, label='LRS2 Standards')
plt.plot_date(np.array(alldT)[inds], percentile_filter(SS, 75, size=50), 'k-', color='darkgray', lw=3, label='Guider')
plt.ylim([0, 1.4])
plt.legend()
plt.xlim([datetime.date(2018, 6, 1), datetime.date(2020, 9, 1)])
plt.gcf().autofmt_xdate()

# Plot formatters
myFmt = mdates.DateFormatter('%m/%d/%y')
weeks = mdates.WeekdayLocator()  # every week
months = mdates.MonthLocator()  # every month
days = mdates.DayLocator()  # every day
mL = MultipleLocator(0.1)
ML = MultipleLocator(0.5)

plt.ylabel('Transparency')
plt.gca().grid('on')
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
plt.gca().tick_params(axis='y', which='both', labelright='on', right=True)
plt.savefig('date_throughput.png', dpi=300)