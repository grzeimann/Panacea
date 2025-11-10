#!/usr/bin/env python3
"""Panacea main script.

What this script does (overview):
- Parses command-line options for date(s), channel(s), and simple switches.
- Gathers calibration frames (bias, twilight/flat, arc) for each IFU channel.
- Derives per-fiber trace, wavelength solution, and flat-field calibrations.
- Loads packaged configuration resources (line lists, DAR tables, responses).
- Locates science exposures and runs a streamlined reduction via
  ``routine.big_reduction`` that performs extraction, sky subtraction,
  optional source finding/extraction, and writes products.
"""
import numpy as np
import os.path as op
import glob
import sys
import argparse as ap
from pathlib import Path
from astropy.io import fits
from scipy.interpolate import interp1d
from importlib import resources

from .input_utils import setup_logging
from .io import (
    get_cal_path,
    get_tarname_from_filename,
    get_filenames_from_tarfolder,
)
from .ccd import (
    get_masterbias,
    get_mastertwi,
    get_bigW,
    get_bigF,
    get_twiflat_field,
    get_masterarc,
)
from .trace import get_trace
from .wavelength import get_wavelength_from_arc
from .routine import get_ifucenfile, big_reduction
from .fiber import get_spectra, weighted_extraction
from .utils import get_objects, check_if_standard, get_config_file, read_arc_lines


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

    parser.add_argument("-cw", "--central_wave",
                        help='''Central Wavelength for collapsed Frame''',
                        type=float, default=None)

    parser.add_argument("-wb", "--wavelength_bin",
                        help='''Wavelength Bin to collapse over (+/- bin size)''',
                        type=float, default=10.)

    parser.add_argument("-re", "--reduce_eng",
                        help='''Reduce Engineer Data''',
                        action="count", default=0)

    parser.add_argument("-rf", "--reduce_flt",
                        help='''Reduce Flat Data''',
                        action="count", default=0)

    parser.add_argument("-rd", "--reduce_drk",
                        help='''Reduce Dark Data''',
                        action="count", default=0)

    parser.add_argument("--baseraw",
                        help='''Base directory containing LRS2 raw data (tarballs). Overrides the built-in default.''',
                        type=str, default='/Users/grz85/data/LRS2')

    args = parser.parse_args(args=None)

    # Sanity check: source_x/source_y require central_wave, and vice versa
    for i in ['source_x', 'source_y']:
        for j in ['source_x', 'source_y', 'central_wave']:
            if i == j:
                continue
            if getattr(args, i) is not None:
                if getattr(args, j) is None:
                    log.error('%s was set but not %s.' % (i, j))
                    sys.exit(1)

    args.sides = [x.replace(' ', '') for x in args.sides.split(',')]

    # Channel metadata: [tag, name, multi-id, wavelength limits, amps, collapse window, arc sequence]
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

    # Select channel descriptors based on requested --sides
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

    # Locate packaged config directory and files
    lrs2config_dir = resources.files('panacea') / 'lrs2_config'

    # Defaults for this environment (paths may be site-specific)
    # Use packaged fplane.txt via importlib.resources and convert to filesystem path
    fplane_file = str(get_config_file('fplane.txt'))
    twi_date = args.date
    sci_date = args.date

    # Instrument tag used in path templates
    instrument = 'lrs2'

    # Placeholder: dither offsets (not currently used in quicklook path)
    dither_pattern = np.zeros((50, 2))

    # Base raw data directory and path patterns for science/calibration discovery
    baseraw = args.baseraw

    # Tarball glob and internal FITS path templates (filled per IFU slot/exp)
    sci_tar = op.join(baseraw, sci_date, '%s', '%s000*.tar')
    sci_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp%s',
                       '%s', '2*_%sLL*sci.fits')
    cmp_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                       '%s', '2*_%sLL_cmp.fits')
    twi_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                       '%s', '2*_%sLL_twi.fits')
    bias_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                        '%s', '2*_%sLL_zro.fits')

    # Optional switches: repoint science pattern to engineering/flat/dark
    if args.reduce_eng:
        sci_path = sci_path.replace('sci', 'eng')

    if args.reduce_flt:
        sci_path = sci_path.replace('sci', 'flt')

    if args.reduce_drk:
        sci_path = sci_path.replace('sci', 'drk')
    # Accumulators (legacy placeholders for extended products/logging)
    fiberpos, fiberspec = ([], [])
    log.info('Beginning the long haul.')
    allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [], [], [])


    # Per-channel processing loop (each entry corresponds to one spectrograph arm)
    for info in listinfo:
        specinit, specname, multi, lims, amps, slims, arc_names = info
        if int(args.date) < 20161101:
            nnn = specname  # '%s_old' % specname
        else:
            nnn = specname
        # Load arc line list via package resources
        with get_config_file(f'lines_{nnn}.dat').open('r') as f:
            arc_lines = read_arc_lines(f)
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
            # Lookup per-amp fiber center positions for spatial mapping
            amppos = get_ifucenfile(specname, amp, lrs2config=str(lrs2config_dir))
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
                                    args.date, lrs2config=str(lrs2config_dir))

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

            log.info('Getting Wavelength for ifuslot, %s, and amp, %s' %
                     (ifuslot, amp))

            wave = get_wavelength_from_arc(masterarc, trace, arc_lines, specname,
                                           amp, int(args.date), otherimage=def_arc)

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
        # Combine per-amp calibrations into single arrays ordered by fiber index
        calinfo = [np.vstack([package[0][i], package[1][i]])
                   for i in np.arange(len(package[0]))]
        masterarc, masterflt, masterFlat = [np.vstack(x) for x in [marc, mtwi, mflt]]
        # Offset trace rows for the upper amp by the number of lower-amp fibers
        calinfo[1][package[0][1].shape[0]:, :] += package[0][2].shape[0]
        log.info('Getting flat for ifuslot, %s, side, %s' % (ifuslot, specname))
        # Build twilight flat field (fiber profile) used for extractions
        twiflat = get_twiflat_field(twifiles, amps, calinfo[0], calinfo[1],
                                    calinfo[2], calinfo[3], specname)
        calinfo.insert(2, twiflat)
        flatspec = get_spectra(calinfo[2], calinfo[1])
        for mfile in [masterarc, masterflt, masterFlat]:
            # Extract per-fiber arc spectra for each calibration image and
            # interpolate them onto the common wavelength grid used downstream.
            masterarcerror = np.sqrt(3. ** 2 + np.where(mfile > 0., mfile, 0.))
            arcspec, ae, Cc, Yyy, Fff = weighted_extraction(
                mfile, masterarcerror, calinfo[2], calinfo[1], cthresh=500
            )
            arc_interp_to_common = np.zeros((calinfo[0].shape[0], len(commonwave)))
            for fiber in np.arange(calinfo[0].shape[0]):
                interp = interp1d(calinfo[0][fiber], arcspec[fiber],
                                  kind='linear', fill_value='extrapolate')
                arc_interp_to_common[fiber] = interp(commonwave)
            calinfo.append(arc_interp_to_common)
        bigF = get_bigF(calinfo[1], calinfo[2])
        calinfo.append(bigF)

        #####################
        # SCIENCE REDUCTION #
        #####################
        # Locate the first exposure (exp01) per science observation for this IFU slot.
        # We then iterate observations and run the per-exposure reduction.
        response = None
        sci_glob = sci_path % (instrument, instrument, '0000*', '01', instrument, ifuslot)
        base_files_nested = []
        for tarname in glob.glob(get_tarname_from_filename(sci_glob)):
            base_files_nested.append(get_filenames_from_tarfolder(tarname, sci_glob))
        flat_list = [item for sublist in base_files_nested for item in sublist]
        base_files = [f for f in sorted(flat_list) if "exp01" in f]

        # Derive 7-digit obsid from folder structure and read object headers
        all_sci_obs = [op.basename(op.dirname(op.dirname(op.dirname(fn))))[-7:]
                       for fn in base_files]
        objects = get_objects(base_files, ['OBJECT', 'EXPTIME'])

        # Load a packaged average response if none provided
        if response is None:
            log.info('Getting average response')
            basename = 'LRS2/CALS'
            with get_config_file(f'response_{specname}.fits').open('rb') as fh:
                R = fits.open(fh)
                response = R[0].data[1] * 1.

        # Assemble a CALS bundle capturing key calibration products for this channel
        hdus = []
        hdu_names = ['wavelength', 'trace', 'flat', 'bigW', 'masterbias',
                     'xypos', 'dead', 'arcspec', 'fltspec', 'Flatspec', 'bigF']
        for i, cal in enumerate(calinfo):
            hdu_class = fits.PrimaryHDU if i == 0 else fits.ImageHDU
            hdus.append(hdu_class(cal))
        hdus.append(fits.ImageHDU(masterarc)); hdu_names.append('masterarc')
        hdus.append(fits.ImageHDU(masterFlat)); hdu_names.append('masterFlat')
        if response is not None:
            hdus.append(fits.ImageHDU(np.array([commonwave, response], dtype=float)))
            hdu_names.append('response')
        for hdu, name in zip(hdus, hdu_names):
            hdu.header['EXTNAME'] = name
        basename = 'LRS2/CALS'
        Path(basename).mkdir(parents=True, exist_ok=True)
        fits.HDUList(hdus).writeto(
            op.join(basename, 'cal_%s_%s.fits' % (args.date, specname)),
            overwrite=True
        )
        # Iterate per observation base file; filter by --object or include standards
        for sci_obs, obj, bf in zip(all_sci_obs, objects, base_files):
            log.info('Checkpoint --- Working on %s, %s' % (bf, specname))
            if args.object is None:
                big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                              ifuslot, specname, response=response, central_wave=args.central_wave,
                              wavelength_bin=args.wavelength_bin, correct_ftf_flag=args.correct_ftf,
                              fplane_file=fplane_file, date_str=args.date)
            else:
                if args.object.lower() in obj[0].lower():
                    big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                                  ifuslot, specname, response=response, central_wave=args.central_wave,
                                  wavelength_bin=args.wavelength_bin, correct_ftf_flag=args.correct_ftf,
                                  fplane_file=fplane_file, date_str=args.date)
                if check_if_standard(obj[0]) and (ifuslot in obj[0]):
                    big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                                  ifuslot, specname, response=response,
                                  central_wave=args.central_wave, wavelength_bin=args.wavelength_bin,
                                  correct_ftf_flag=args.correct_ftf,
                                  fplane_file=fplane_file, date_str=args.date)


if __name__ == '__main__':
    main()