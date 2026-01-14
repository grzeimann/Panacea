# Configuration and Packaged Resources

Panacea bundles instrument configuration resources with the Python package under panacea/lrs2_config. At runtime, code locates these files via importlib.resources, so you normally do not need to manage absolute paths.

What’s inside lrs2_config repository

1) IFU fiber-center mappings (used for spatial geometry)
- LRS2_B_UV_mapping.txt
- LRS2_B_OR_mapping.txt
- LRS2_R_NR_mapping.txt
- LRS2_R_FR_mapping.txt

Where used: routine.get_ifucenfile(side, amp) loads these and returns (x, y) fiber coordinates for each amplifier half.

2) DAR reference tables (used to model differential atmospheric refraction)
- dar_BL.dat
- dar_BR.dat
- dar_RL.dat
- dar_RR.dat

Where used: routine.big_reduction reads the channel-appropriate table (BL for uv, BR for orange, RL for red, RR for farred) to build wavelength-dependent DAR offsets (xoff, yoff).

3) Arc line lists (used for wavelength calibration)
- lines_uv.dat
- lines_orange.dat
- lines_red.dat
- lines_farred.dat

Where used: run_panacea.py opens lines_{channel}.dat and passes them to wavelength.get_wavelength_from_arc, which detects arc peaks and fits dispersion per fiber.

4) Average response functions (packaged medians; used if no user-provided response)
- response_uv.fits
- response_orange.fits
- response_red.fits
- response_farred.fits

Where used: run_panacea.py loads response_{channel}.fits when building the per-channel CALS bundle and to optionally scale spectra. Note: There are also response_BL.dat, response_BR.dat, response_RL.dat, response_RR.dat text files present; these are legacy and not currently read by the quicklook path.

5) Skyline masks (used to suppress skylines during source finding and statistics)
- orange_skylines.dat
- red_skylines.dat
- farred_skylines.dat
- uv channel: there is no uv_skylines.dat in this tree. The code tolerates a missing file and applies a small built-in mask for uv in utils.mask_skylines_cosmics.

Where used: utils.mask_skylines_cosmics loads f"{channel}_skylines.dat" if present.

6) Focal-plane geometry (used for astrometric mapping)
- fplane.txt

Where used: Astrometry/FPlane classes in astrometry.py. The quicklook reduction passes this into routine.big_reduction so reconstructed cubes can be mapped into IFU plane coordinates and (if headers allow) sky coordinates.

7) Fiber_Locations (reference fiber traces; pick closest date)
- Directory structure: Fiber_Locations/YYYYMMDD/fiber_loc_{specid}_{ifuslot}_{ifuid}_{amp}.txt
- Present dates: 20160724, 20170305, 20181108 with files for amps LL, LU, RL, RU and specids 501/502/503.

Where used: trace.get_trace_reference finds the closest-in-time subdirectory and loads the matching file; trace.get_trace uses it to seed and stabilize trace finding.

8) Additional calibration aids currently packaged but not used by run/routine
- ftf_BL.fits, ftf_BR.fits, ftf_RL.fits, ftf_RR.fits (candidate fiber-to-fiber templates)
- uv_wavelength.fits, orange_wavelength.fits, red_wavelength.fits, farred_wavelength.fits (candidate wavelength priors)

Overriding the packaged defaults
- Within Python: You can open and use your own files, but for the CLI quicklook, the simplest “override” is to place a file with the same name in the installed package data. This is easiest in an editable/development install. Example for replacing the ORANGE arc lines:
  cp my_lines_orange.dat src/panacea/lrs2_config/lines_orange.dat
  (Reinstall or rebuild the package as needed.)

File formats (summary)
- Mapping files (LRS2_*_mapping.txt): ASCII with columns [fiber_id, x, y, ...]; routine.get_ifucenfile reads columns 0–2 (id, x, y) and splits by amplifier.
- DAR tables (dar_*.dat): ASCII with header; columns wave, x_0, y_0 (in IFU pixel units) sampled across the channel bandpass.
- Arc line lists (lines_*.dat): ASCII, whitespace-separated with optional comments (#). utils.read_arc_lines parses the first four columns as [wavelength, approx_x, relative_intensity, name].
- Response FITS (response_*.fits): Primary HDU contains a 2-row array [wavelength, response] used to scale spectra.
- Skyline masks (*_skylines.dat): One wavelength per line; utils.mask_skylines_cosmics masks ±~6 Å around each listed wavelength.
- fplane.txt: Whitespace-separated columns describing IFU slot, positions, spectrograph id/slot, ifuid, rotation, plate scale. astrometry.FPlane expects 8 fields per row and skips comments (#).
- Fiber_Locations files: Per-amp reference arrays with per-fiber [column_index, dead_flag]; trace.get_trace_reference reads them as text.

See also
- CLI Usage: ../user-guide/cli.md
- API Reference for utils.get_config_file and related helpers: ../api/index.md
