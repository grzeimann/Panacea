#!/usr/bin/env python3
"""Panacea CLI wrapper.

This script delegates to the package entry point. It allows direct execution
from a source checkout without installation: `python scripts/run_panacea.py -h`.
"""
import numpy as np
import os.path as op
import glob
import sys
import argparse as ap
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d
from input_utils import setup_logging

# Import migrated functions from the Panacea package
from panacea.io import (
    get_cal_path,
    get_tarname_from_filename,
    get_filenames_from_tarfolder,
)
from panacea.ccd import (
    get_masterbias,
    get_mastertwi,
    get_bigW,
    get_bigF,
    get_twiflat_field,
    get_masterarc,
)
from panacea.trace import get_trace
from panacea.wavelength import get_wavelength_from_arc
from panacea.routine import get_ifucenfile, big_reduction
from panacea.fiber import get_spectra, weighted_extraction
from panacea.utils import get_script_path, get_objects, check_if_standard


standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4',
                  'BD+40_4032', 'HILTNER_102',
                  'BD_+26_2606', 'GD_248', 'FEIGE_56', 'FEIGE_92',
                  'HZ_15', 'FEIGE_98', 'BD+08_2015', 'BD+25_3941',
                  'FEIGE_15', 'FEIGE_25', 'SA_95-42', 'BD+28_4211',
                  'HR6203']

log = setup_logging('panacea_quicklook')


def main():
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date for reduction''',
                        type=str, default='20181108')

    parser.add_argument("-s", "--sides",
                        help='''"uv,orange,red,farred"''',
                        type=str, default="uv,orange,red,farred")

    parser.add_argument("-o", "--object",
                        help='''Object name, no input reduces all objects''',
                        type=str, default=None)

    parser.add_argument("-uf", "--use_flat",
                        help='''Use FLT instead of Twi''',
                        action="count", default=0)

    parser.add_argument("-cf", "--correct_ftf",
                        help='''Correct fiber to fiber''',
                        action="count", default=0)

    parser.add_argument("-md", "--model_dar",
                        help='''model DAR''',
                        action="count", default=0)

    parser.add_argument("-cw", "--central_wave",
                        help='''Central Wavelength for collapsed Frame''',
                        type=float, default=None)

    parser.add_argument("-wb", "--wavelength_bin",
                        help='''Wavelength Bin to collapse over (+/- bin size)''',
                        type=float, default=10.)

    parser.add_argument("-sx", "--source_x",
                        help='''Source's x position at the central_wave''',
                        type=float, default=None)

    parser.add_argument("-sy", "--source_y",
                        help='''Source's y position at the central_wave''',
                        type=float, default=None)

    parser.add_argument("-ssd", "--standard_star_date",
                        help='''Standard Star Date for response function,
                        example: 20181101''',
                        type=str, default=None)

    parser.add_argument("-sso", "--standard_star_obsid",
                        help='''Standard Star ObsID for response function,
                        example: 0000012''',
                        type=str, default=None)

    parser.add_argument("-ad", "--arc_date",
                        help='''Arc Date for reduction''',
                        type=str, default=None)

    parser.add_argument("-td", "--twi_date",
                        help='''Twilight Date for reduction''',
                        type=str, default=None)

    parser.add_argument("-re", "--reduce_eng",
                        help='''Reduce Engineer Data''',
                        action="count", default=0)

    parser.add_argument("-rf", "--reduce_flt",
                        help='''Reduce Flat Data''',
                        action="count", default=0)

    parser.add_argument("-rd", "--reduce_drk",
                        help='''Reduce Dark Data''',
                        action="count", default=0)

    args = parser.parse_args(args=None)

    if args.standard_star_obsid is not None:
        args.standard_star_obsid = '%07d' % int(args.standard_star_obsid)
        if args.standard_star_date is None:
            log.error('Please include --standard_star_date DATE with call.')

    for i in ['source_x', 'source_y']:
        for j in ['source_x', 'source_y', 'central_wave']:
            if i == j:
                continue
            if getattr(args, i) is not None:
                if getattr(args, j) is None:
                    log.error('%s was set but not %s.' % (i, j))
                    sys.exit(1)

    args.sides = [x.replace(' ', '') for x in args.sides.split(',')]

    blueinfo = [['BL', 'uv', '503_056_7001', [3640., 4645.], ['LL', 'LU'],
                 [4350., 4375.], ['hg_b', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']],
                ['BR', 'orange', '503_056_7001',
                 [4635., 6950.], ['RU', 'RL'], [6270., 6470.],
                 ['hg_b', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']]]
    redinfo = [['RL', 'red', '502_066_7002', [6450., 8400.], ['LL', 'LU'],
                [7225., 7425.], ['hg_r', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']],
               ['RR', 'farred', '502_066_7002',
                [8275., 10500.], ['RU', 'RL'], [9280., 9530.],
                ['hg_r', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']]]

    listinfo = []
    for side in args.sides:
        if side.lower() == 'uv':
            listinfo.append(blueinfo[0])
        if side.lower() == 'orange':
            listinfo.append(blueinfo[1])
        if side.lower() == 'red':
            listinfo.append(redinfo[0])
        if side.lower() == 'farred':
            listinfo.append(redinfo[1])

    fplane_file = '/work/03730/gregz/maverick/fplane.txt'
    twi_date = args.date
    sci_date = args.date

    # FOR LRS2
    instrument = 'lrs2'

    dither_pattern = np.zeros((50, 2))

    baseraw = '/work/03946/hetdex/maverick'

    sci_tar = op.join(baseraw, sci_date, '%s', '%s000*.tar')
    sci_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp%s',
                       '%s', '2*_%sLL*sci.fits')
    cmp_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                       '%s', '2*_%sLL_cmp.fits')
    twi_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                       '%s', '2*_%sLL_twi.fits')
    bias_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                        '%s', '2*_%sLL_zro.fits')

    if args.reduce_eng:
        sci_path = sci_path.replace('sci', 'eng')

    if args.reduce_flt:
        sci_path = sci_path.replace('sci', 'flt')

    if args.reduce_drk:
        sci_path = sci_path.replace('sci', 'drk')
    # LRS2-R
    fiberpos, fiberspec = ([], [])
    log.info('Beginning the long haul.')
    allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [],
                                                               [], [])

    DIRNAME = get_script_path()

    for info in listinfo:
        specinit, specname, multi, lims, amps, slims, arc_names = info
        if int(args.date) < 20161101:
            nnn = specname  # '%s_old' % specname
        else:
            nnn = specname
        arc_lines = Table.read(op.join(DIRNAME, 'lrs2_config/lines_%s.dat' %
                                       nnn), format='ascii')
        commonwave = np.linspace(lims[0], lims[1], 2064)
        specid, ifuslot, ifuid = multi.split('_')
        package = []
        marc, mtwi, mflt = ([], [], [])
        twipath = twi_path % (instrument, instrument, '00000*', instrument,
                              ifuslot)
        twifiles = get_cal_path(twipath, args.date, ndays=15)
        flt_path = (twi_path.replace('twi', 'flt') %
                    (instrument, instrument, '00000*', instrument, ifuslot))
        fltfiles = get_cal_path(flt_path, args.date, ndays=15)
        if args.use_flat:
            twifiles = fltfiles
        for amp in amps:
            amppos = get_ifucenfile(specname, amp)
            ##############
            # MASTERBIAS #
            ##############
            log.info('Getting Masterbias for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))
            zro_path = bias_path % (instrument, instrument, '00000*', instrument,
                                    ifuslot)
            zrofiles = get_cal_path(zro_path, args.date, ndays=2)
            masterbias = get_masterbias(zrofiles, amp)

            #####################
            # MASTERTWI [TRACE] #
            #####################
            log.info('Getting MasterFlat for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))
            log.info('Getting Trace for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))
            log.info('Number of twi files: %i' % len(twifiles))

            masterflt = get_mastertwi(twifiles, amp, masterbias)
            trace, dead = get_trace(masterflt, specid, ifuslot, ifuid, amp,
                                    args.date)

            ##########################
            # MASTERARC [WAVELENGTH] #
            ##########################
            log.info('Getting MasterArc for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))
            lamp_path = cmp_path % (instrument, instrument, '00000*', instrument,
                                    ifuslot)
            lampfiles = get_cal_path(lamp_path, args.date, ndays=15)
            log.info('Number of arc files: %i' % len(lampfiles))
            masterarc = get_masterarc(lampfiles, amp, arc_names, masterbias,
                                      specname, trace)

            lampfiles = get_cal_path(lamp_path.replace(args.date, '20181201'),
                                     '20181201', ndays=3)
            def_arc = get_masterarc(lampfiles, amp,
                                    arc_names, masterbias, specname, trace)

            # fits.PrimaryHDU(masterarc).writeto('/work/03946/hetdex/maverick/run_lrs2/wtf_%s_%s.fits' % (ifuslot, amp), overwrite=True)
            log.info('Getting Wavelength for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))

            wave = get_wavelength_from_arc(masterarc, trace, arc_lines, specname,
                                           amp, int(args.date), otherimage=def_arc)
            # fits.PrimaryHDU(wave).writeto('test_wave.fits', overwrite=True)

            #################################
            # TWILIGHT FLAT [FIBER PROFILE] #
            #################################
            log.info('Getting bigW for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))
            bigW = get_bigW(amp, wave, trace, masterbias)
            package.append([wave, trace, bigW, masterbias, amppos, dead])
            log.info('Number of flt files: %i' % len(fltfiles))

            masterFlat = get_mastertwi(fltfiles, amp, masterbias)

            marc.append(masterarc)
            mtwi.append(masterflt)
            mflt.append(masterFlat)
        # Normalize the two amps and correct the flat
        calinfo = [np.vstack([package[0][i], package[1][i]])
                   for i in np.arange(len(package[0]))]
        masterarc, masterflt, masterFlat = [np.vstack(x) for x in [marc, mtwi, mflt]]
        calinfo[1][package[0][1].shape[0]:, :] += package[0][2].shape[0]
        log.info('Getting flat for ifuslot, %s, side, %s' % (ifuslot, specname))
        twiflat = get_twiflat_field(twifiles, amps, calinfo[0], calinfo[1],
                                    calinfo[2], calinfo[3], specname)
        calinfo.insert(2, twiflat)
        flatspec = get_spectra(calinfo[2], calinfo[1])
        for mfile in [masterarc, masterflt, masterFlat]:
            masterarcerror = np.sqrt(3. ** 2 + np.where(mfile > 0., mfile, 0.))
            arcspec, ae, Cc, Yyy, Fff = weighted_extraction(mfile, masterarcerror,
                                                            calinfo[2], calinfo[1],
                                                            cthresh=500)
            sP = np.zeros((calinfo[0].shape[0], len(commonwave)))
            for fiber in np.arange(calinfo[0].shape[0]):
                I = interp1d(calinfo[0][fiber], arcspec[fiber],
                             kind='linear', fill_value='extrapolate')
                sP[fiber] = I(commonwave)
            calinfo.append(sP)
        bigF = get_bigF(calinfo[1], calinfo[2])
        calinfo.append(bigF)
        #####################
        # SCIENCE REDUCTION #
        #####################
        response = None
        pathS = sci_path % (instrument, instrument, '0000*',
                            '01', instrument, ifuslot)
        basefiles = []
        for tarname in glob.glob(get_tarname_from_filename(pathS)):
            basefiles.append(get_filenames_from_tarfolder(tarname, pathS))
        flat_list = [item for sublist in basefiles for item in sublist]
        basefiles = [f for f in sorted(flat_list) if "exp01" in f]

        all_sci_obs = [op.basename(op.dirname(op.dirname(op.dirname(fn))))[-7:]
                       for fn in basefiles]
        objects = get_objects(basefiles, ['OBJECT', 'EXPTIME'])
        if response is None:
            log.info('Getting average response')
            basename = 'LRS2/CALS'
            R = fits.open(op.join(DIRNAME,
                                  'lrs2_config/response_%s.fits' % specname))
            response = R[0].data[1] * 1.

        f = []
        names = ['wavelength', 'trace', 'flat', 'bigW', 'masterbias',
                 'xypos', 'dead', 'arcspec', 'fltspec', 'Flatspec', 'bigF']
        for i, cal in enumerate(calinfo):
            if i == 0:
                func = fits.PrimaryHDU
            else:
                func = fits.ImageHDU
            f.append(func(cal))
        f.append(fits.ImageHDU(masterarc))
        names.append('masterarc')
        f.append(fits.ImageHDU(masterFlat))
        names.append('masterFlat')
        if response is not None:
            f.append(fits.ImageHDU(np.array([commonwave, response], dtype=float)))
            names.append('response')
        for fi, n in zip(f, names):
            fi.header['EXTNAME'] = n
        basename = 'LRS2/CALS'
        Path(basename).mkdir(parents=True, exist_ok=True)

        fits.HDUList(f).writeto(op.join(basename,
                                        'cal_%s_%s.fits' % (args.date, specname)),
                                overwrite=True)
        for sci_obs, obj, bf in zip(all_sci_obs, objects, basefiles):
            log.info('Checkpoint --- Working on %s, %s' % (bf, specname))
            if args.object is None:
                big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                              ifuslot, specname, response=response,
                              central_wave=args.central_wave, wavelength_bin=args.wavelength_bin,
                              source_x=args.source_x, source_y=args.source_y,
                              correct_ftf_flag=bool(args.correct_ftf))
            else:
                if args.object.lower() in obj[0].lower():
                    big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                                  ifuslot, specname, response=response,
                                  central_wave=args.central_wave, wavelength_bin=args.wavelength_bin,
                                  source_x=args.source_x, source_y=args.source_y,
                                  correct_ftf_flag=bool(args.correct_ftf))
                if check_if_standard(obj[0]) and (ifuslot in obj[0]):
                    big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                                  ifuslot, specname, response=response,
                                  central_wave=args.central_wave, wavelength_bin=args.wavelength_bin,
                                  source_x=args.source_x, source_y=args.source_y,
                                  correct_ftf_flag=bool(args.correct_ftf))

if __name__ == "__main__":
    main()
