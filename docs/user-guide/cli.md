# CLI Usage

This page documents the panacea-lrs2 command-line interface.

Run help:
```bash
panacea-lrs2 -h
```

Common options (as implemented in the current CLI):
```bash
usage: panacea-lrs2 [-h] [-d DATE] [-s SIDES] [-o OBJECT] [-uf] [-cf] [-cw CENTRAL_WAVE] [-wb WAVELENGTH_BIN] [-re] [-rf] [-rd] [--baseraw BASERAW]usage: panacea-lrs2 [-h] [-d DATE] [-s SIDES] [-o OBJECT] [-uf] [-cf] [-cw CENTRAL_WAVE] [-wb WAVELENGTH_BIN] [-re] [-rf] [-rd] [--baseraw BASERAW]

options:
  -h, --help            show this help message and exit
  -d, --date DATE       Observation date YYYYMMDD
  -s, --sides SIDES     Comma-separated channels or just a single channel
  -o, --object OBJECT   Substring to match OBJECT header (omit to reduce all)
  -uf, --use_flat       Use internal flat (FLT) instead of twilight flats
  -cf, --correct_ftf    Enable additional fiber-to-fiber correction using sky emission
  -cw, --central_wave CENTRAL_WAVE
                        Center wavelength (Å) for collapsed image
  -wb, --wavelength_bin WAVELENGTH_BIN
                        Half-width (Å) of collapse window around center
  -re, --reduce_eng     Use engineering (ENG) exposures instead of SCI
  -rf, --reduce_flt     Use flat (FLT) frames as science
  -rd, --reduce_drk     Use dark (DRK) frames as science
  --baseraw BASERAW     Base directory containing LRS2 raw data (tarballs). Overrides the built-in default.
```

Notes:
- The collapse window is [central_wave - wavelength_bin, central_wave + wavelength_bin]. If central_wave is not provided, the center of the channel’s wavelength range is used.
- If you provide central_wave, do not try to reduce multiple channels as that should be done for one channel at a time
- Flags using action "count" (-uf, -cf, -re, -rf, -rd) are treated as booleans: include the flag to enable the behavior.

Examples:
- Reduce only the orange channel for a given night (suggested usage):
```bash
panacea-lrs2 -d 20181108 -s orange
```
- Reduce all channels and enable ftf correction:
```bash
panacea-lrs2 -d 20181108  -cf
```
- Reduce all four channels using internal flats 
```bash
panacea-lrs2 -d 20181108 -uf
```

See also: [Quickstart](../getting-started/quickstart.md), [Running on TACC](../tacc/running.md)
