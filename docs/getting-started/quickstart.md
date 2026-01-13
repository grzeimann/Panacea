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
panacea-lrs2 -d 20251016 -s farred

# Reduce UV for a specific date
panacea-lrs2 -d 20181108 -s uv
```

## Notes for local runs
- The quicklook runner expects HET raw data in a TACC-like directory structure (tarballs with standard internal paths).
- Use the --baseraw option to point Panacea to your local mirror of the LRS2 raw data.
  Example: `--baseraw /data/LRS2`
- Packaged configuration files (line lists, DAR tables, fplane.txt, responses) are bundled with the package and found automatically via importlib.resources.

See also: [CLI Usage](../user-guide/cli.md) and [Configuration](../user-guide/configuration.md)
