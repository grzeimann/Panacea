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
from telluricabs import TelluricAbs
import os.path as op
from input_utils import setup_basic_parser, setup_logging


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

args = parser.parse_args(['--instr', 'lrs2', '--rootdir',
                          '/Users/gregz/cure/reductions', '--side', 'red',
                          '-d', '20171126', '-o', '0000022', '-e', '1'])
args.log = setup_logging(logname='response_lrs2')
args.observation = int(args.observation)
args.exposure_number = int(args.exposure_number)

fplane_file = ('/work/03946/hetdex/maverick/virus_config/fplane/'
               'fplane20180419.txt')

fig9, ax9 = plt.subplots(1, 1, figsize=(8, 6))
if args.side == 'blue':
    multiname = 'multi_503_056_7001'
    sides = ['BL', 'BR']
else:
    multiname = 'multi_502_066_7002'
    sides = ['RL', 'RR']

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
#    T = TelluricAbs(P.dar.rect_wave, P.dar.flux, 20., 300., 797., 22.)
#    absorp = T.get_model(P.dar.rect_wave.min()-5., P.dar.rect_wave.max()+5.,
#                         P.dar.rect_wave, resolution=1.05, o2=210000.,
#                         humidity=30)

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
#    P.clam = P.clam / absorp
    Alam = P.get_extinction_mag(np.mean(P.ra), np.mean(P.dec),
                                P.dar.rect_wave / 1e4)
    P.clam = 10**(0.4 * Alam) * P.clam
    P.compare_spectrum_to_standard()
    ax9.plot(P.R_wave, P.smooth_R)
    ax9.set_ylim([0, 0.6e-8])
    ax9.set_xlim([3000, 11000])

fig9.savefig('Response_%s_%07d_%s.png' % (args.date, args.observation,
                                          args.side))
