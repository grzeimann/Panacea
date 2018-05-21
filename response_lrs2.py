# -*- coding: utf-8 -*-
"""
Measure response LRS2

@author: gregz
"""
import matplotlib
matplotlib.use('agg')
from reducelrs2 import ReduceLRS2
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import os.path as op
from astropy.io import fits
from scipy.signal import savgol_filter
from input_utils import setup_basic_parser, setup_logging
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import Pipeline


standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B', 'FEIGE_25',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4']


def fit_continuum(wv, sky, skip=3, fil_len=95, func=np.array):
    skym_s = 1. * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    for i in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        outlier = func(sky - sky_sm) > 1.5 * mad
        sel = np.where(outlier)[0]
        for j in np.arange(1, skip+1):
            sel = np.union1d(sel, sel + 1)
            sel = np.union1d(sel, sel - 1)
        sel = np.sort(np.unique(sel))
        sel = sel[skip:-skip]
        good = np.setdiff1d(allind, sel)
        skym_s = 1.*sky
        skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    return sky_sm


def build_filename(args, multiname):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    filename = op.join(args.rootdir, args.date, args.instrument,
                       '%s%07d' % (args.instrument, args.observation),
                       'exp%02d' % args.exposure_number, args.instrument,
                       multiname)
    return filename

parser = setup_basic_parser()

# ['-d', '20171126', '-o', '0000021', '-e', '1',  '-r', '/Users/gregz/cure/reductions', '--side', 'blue', '--instr', 'lrs2']
args = parser.parse_args()
args.log = setup_logging(logname='response_lrs2')
args.observation = int(args.observation)
args.exposure_number = int(args.exposure_number)

fplane_file = ('/work/03946/hetdex/maverick/virus_config/fplane/'
               'fplane20180419.txt')

if args.side == 'blue':
    multiname = 'multi_503_056_7001'
    sides = ['BL', 'BR']
    skipv = 3
else:
    multiname = 'multi_502_066_7002'
    sides = ['RL', 'RR']
    skipv = 7

for side in sides:
    filebase = build_filename(args, multiname)
    outfolder = op.dirname(filebase)
    P = ReduceLRS2(filebase, side, fplane_file=fplane_file)
    for standard in standard_names:
        if standard.lower() in P.object.lower():
            P.object = standard.lower()
    P.dar.rectify(minwave=P.wave_lims[0], maxwave=P.wave_lims[1])
    P.dar.measure_dar()
    P.dar.psfextract()
    for i, ind in enumerate(np.arange(100, len(P.dar.rect_wave), 100)):
        outname = op.join(outfolder, 'psf_standard_%04d_%s.png' % (ind, side))
        P.dar.check_psf_fit(index=ind, outname=outname)

    attr = P.dar.tinker_params + ['fwhm']
    fig, ax = plt.subplots(len(attr), 1, sharex=True, figsize=(len(attr),
                           len(attr)*3))
    T = np.zeros((len(P.dar.dar_wave), len(attr)+1))
    T[:, 0] = P.dar.dar_wave
    for a, at, cnt in zip(ax, attr, np.arange(1, len(attr)+1)):
        a.scatter(P.dar.dar_wave, P.dar.A[:, cnt])
        T[:, cnt] = P.dar.A[:, cnt]
        a.plot(P.dar.dar_wave, getattr(P.dar, 'dar_'+at), color='darkorange')
        a.set_title(at)
    T = Table(T)
    T.write(op.join(outfolder, 'dar_info_%s.dat' % side), format='ascii')
    fig.savefig(op.join(outfolder,
                        'dar_standard_%s_%07d_%s.png' % (args.date,
                                                         args.observation,
                                                         side)))
    plt.close(fig)
    P.convert_units()
    P.clam_old = 1. * P.clam
    x = P.dar.rect_wave
    P.clam = -1. * fit_continuum(P.dar.rect_wave, -P.clam, skip=skipv)
    P.get_standard_spectrum_from_file()
    xl = np.searchsorted(P.standard_wave, P.dar.rect_wave.min())-3
    xh = np.searchsorted(P.standard_wave, P.dar.rect_wave.max())+3
    x1 = P.standard_wave[xl:xh]
    y1 = P.standard_flam[xl:xh]
    y2 = P.dar.robust_poly_fit(x1, y1, 3)
    y3 = np.interp(P.dar.rect_wave, x1, y2)
    area = P.area / 55e4 * (np.pi * 500**2)
    args.log.info('Exposure time: %0.2f' % P.exptime)
    exp = y1 * area * (P.exptime-3.5) / 6.63e-27 / (3e18 / x1)
    x0 = P.dar.rect_wave
    P.flux_binned_wave = x1[1:-1]
    P.flux_bin = np.zeros(x1[1:-1].shape)
    for i in np.arange(1, len(x1)-1):
        low = x1[i-1] / 2. + x1[i] / 2.
        high = x1[i] / 2. + x1[i+1] / 2.
        sel = np.where((x0 > low) * (x0 <= high))[0]
        P.flux_bin[i-1] = P.dar.flux[sel].mean()
    X = P.flux_bin / exp[1:-1]
    sel = np.isfinite(X)
    x = P.flux_binned_wave[sel]
    y = X[sel]
    model = Pipeline([('poly', PolynomialFeatures(degree=7)),
                      ('linear', RANSACRegressor())])
    model.fit(x[:, np.newaxis], y)
    roughthrough = np.interp(P.dar.rect_wave, x, y)
    through = model.predict(P.dar.rect_wave[:, np.newaxis])
    P.R = np.interp(P.dar.rect_wave, x1, y2) / P.clam
    hdu = fits.PrimaryHDU(np.array([P.dar.rect_wave, P.R, P.clam_old,
                                    roughthrough, through], dtype='float32'))
    for key in P.header:
        hdu.header[key] = P.header[key]
    hdu.writeto(op.join(outfolder, 'responsecurve_%s.fits' % side),
                overwrite=True)
