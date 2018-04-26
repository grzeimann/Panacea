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
    P.get_dar_model()
    P.dar.measure_dar(fixed_list=['alpha', 'gamma', 'ratio'])
    P.dar.measure_dar(fixed_list=['x_0', 'y_0'])
    for i, ind in enumerate(np.arange(100, len(P.dar.rect_wave), 100)):
        outname = op.join(outfolder, 'psf_standard_%04d_%s.png' % (ind, side))
        P.dar.check_psf_fit(index=ind, outname=outname)

    attr = P.dar.tinker_params + ['fwhm']
    fig, ax = plt.subplots(len(attr), 1, sharex=True, figsize=(len(attr),
                           len(attr)*3))
    cnt = 0
    for a, at in zip(ax, attr):
        a.scatter(P.dar.dar_wave, P.dar.A[:, cnt+1])
        a.plot(P.dar.dar_wave, getattr(P.dar, 'dar_'+at), color='darkorange')
        a.set_title(at)
        cnt = cnt+1
    fig.savefig(op.join(outfolder,
                        'dar_standard_%s_%07d_%s.png' % (args.date,
                                                         args.observation,
                                                         side)))
    plt.close(fig)
    P.convert_units()
    P.clam_old = 1. * P.clam
    x = P.dar.rect_wave
    P.clam_unred = -1. * fit_continuum(P.dar.rect_wave, -P.clam, skip=skipv)
    Alam = P.get_extinction_mag(np.mean(P.ra), np.mean(P.dec),
                                P.dar.rect_wave / 1e4)

    P.clam = 10**(0.4 * Alam) * P.clam_unred
    P.get_standard_spectrum_from_file()
    xl = np.searchsorted(P.standard_wave, P.dar.rect_wave.min())-2
    xh = np.searchsorted(P.standard_wave, P.dar.rect_wave.max())+2
    x1 = P.standard_wave[xl:xh]
    y1 = P.standard_flam[xl:xh]
    y2 = P.dar.robust_poly_fit(x1, y1, 3)
    P.R = np.interp(P.dar.rect_wave, x1, y2) / P.clam
    hdu = fits.PrimaryHDU(np.array([P.dar.rect_wave, P.R, P.clam_old,
                                    P.clam_unred], dtype='float32'))
    hdu[0].header = P.header
    hdu.writeto(op.join(outfolder, 'responsecurve_%s.fits' % side),
                overwrite=True)
