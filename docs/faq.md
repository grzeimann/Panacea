# Frequently Asked Questions

This page collects common questions about Panacea usage, calibration, and outputs.

## What wavelength units do Panacea products use?
- Wavelengths are in vacuum Angstroms (Å).
- Spectra are rectified onto a common vacuum wavelength grid per channel during extraction.

## How is the flux calibration handled? Do I need a nightly standard star?
- By default, Panacea applies a packaged average response curve for each LRS2 channel (uv, orange, red, farred). It does not require or use a standard star from the same night.
- Why this works: The Hobby–Eberly Telescope (HET) operates at a fixed altitude and typically observes at a near-constant airmass. This stability makes a median “standard response” effective for quicklook calibration.
- Practical implications:
  - The response places spectra on a consistent relative flux scale suitable for quicklook science and comparisons across nights.
  - The initial calibration is good to 10% across wavelengths in relative calibration and ~20% in absolute.
  - If you need per-night absolute spectrophotometry, you can derive a bespoke response outside of the quicklook path and apply it post‑hoc.
  - LRS2Multi (below) is great for secondary calibration

## Are the outputs single-channel or combined across channels?
- Panacea’s reduction products are per channel (one of uv, orange, red, farred) and are written separately.
- To stitch/combine channels into a single spectrum or cube, use LRS2Multi:
  - LRS2Multi repository: https://github.com/grzeimann/LRS2Multi

## Can I select which wavelengths are emphasized in the collapsed image?
- Yes. Use the CLI options:
  - --central_wave: Center wavelength (Å) for the collapsed image window.
  - --wavelength_bin: Half-width (Å) around the center used for collapsing.
- See the CLI docs for details and examples.

## Where can I learn more about configuration files and defaults?
- See the Configuration guide, which lists the packaged line lists, DAR tables, response curves, and geometry resources.

See also:
- [CLI Usage](./user-guide/cli.md)
- [Configuration](./user-guide/configuration.md)
- [Data Products](./data-products/overview.md)
- [Installation](./getting-started/installation.md)
