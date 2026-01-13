# Running on TACC

This page describes interactive and batch ways to run Panacea on TACC.

## Interactive exploration
Start a development session and inspect CLI help:
```bash
idev
panacea-lrs2 -h
```

Common options:
- `-d, --date YYYYMMDD` date of observations (required)
- `-s, --sides uv,orange,red,farred` which channels to reduce (comma-separated)
- `-o, --object NAME` match target OBJECT header containing NAME
- `--use_flat`, `--correct_ftf`
- `--central_wave`, `--wavelength_bin`, `--source_x`, `--source_y`
- `--standard_star_date`, `--standard_star_obsid`

Example:
```bash
panacea-lrs2 -d 20181108 -s uv -o HD_19445
```

## Batch submission
To run reductions in batch for a specific target and date on all four channels:
```bash
cdw
cp /work/03946/hetdex/maverick/run_lrs2/runlrs2general .
runlrs2general DATE TARGET_NAME
```

You will receive a job number upon submission. Monitor logs such as `reductionlrs2daily.oXXXXXX` for progress.

See also: [TACC Overview](./overview.md), [CLI Usage](../user-guide/cli.md)