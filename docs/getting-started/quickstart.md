# Quickstart

This page shows first commands to verify your installation and run basic reductions locally. For running on TACC, see ../tacc/running.md.

## Verify installation
```bash
panacea-lrs2 -h
```

## Example commands
Adjust arguments to your data and environment.
```bash
# Reduce Far Red for a specific date
panacea-lrs2 -d 20260115 -s farred

# Reduce UV for a specific date
panacea-lrs2 -d 20260115 -s uv
```

## Run on sample data (~400 MB)
If new users need to try a real reduction quickly, a small sample dataset is available. Replace DATA_FOLDER and OUTPUT_FOLDER with your paths.

```bash
# 1) Download the sample tarball (about 400 MB)
# Browse to: https://web.corral.tacc.utexas.edu/hetdex/LRS2_test_data/
# Download file: lrs2_20260115_test.tar.gz

# 2) Move it into your chosen data directory
mv lrs2_20260115_test.tar.gz DATA_FOLDER/.

# 3) Extract (and remove the tarball to save space)
cd DATA_FOLDER
tar -xvzf lrs2_20260115_test.tar.gz && rm lrs2_20260115_test.tar.gz

# 4) Run the pipeline from your desired output location (or any directory in the same environment)
cd OUTPUT_FOLDER
# Example: reduce the Orange channel for this date using your data folder as --baseraw
panacea-lrs2 -d 20260115 --baseraw DATA_FOLDER -s orange
```

What to expect
- Output directory LRS2/PROGRAM-ID will be created under your current working directory.
- For this sample, PROGRAM-ID values will be CALS and ORPHANS.
- See [Data products overview](../data-products/overview.md) for details on files and structure.

## Notes for local runs
- The quicklook runner expects HET raw data in a TACC-like directory structure (tarballs with standard internal paths).
- Use the --baseraw option to point Panacea to your local mirror of the LRS2 raw data.
  Example: `--baseraw /data/LRS2`
- Packaged configuration files (line lists, DAR tables, fplane.txt, responses) are bundled with the package and found automatically via importlib.resources.

See also: [CLI Usage](../user-guide/cli.md), [Configuration](../user-guide/configuration.md), and [Tests overview](../development/tests-overview.md)
