---
title: 'Panacea: The LRS2 Data Reduction Pipeline for the Hobby–Eberly Telescope'
tags:
  - astronomy
  - spectroscopy
  - data reduction
  - IFU
  - Hobby–Eberly Telescope
authors:
  - name: Greg Zeimann
    orcid: 0000-0003-2307-0629
    affiliation: 1
affiliations:
  - name: Hobby-Eberly Telescope, University of Texas, Austin
    index: 1
date: 2025-10-20
bibliography: paper.bib
---

# Summary

**Panacea** is the production data-reduction pipeline for the *Low Resolution Spectrograph 2 (LRS2)* on the **Hobby–Eberly Telescope (HET)**. LRS2 is a fiber-fed integral-field spectrograph mounted on the 10-m HET, designed to capture spatially resolved optical spectra across a broad wavelength range. 

The instrument’s four spectrograph channels (UV, Orange, Red, and Far Red) together span 3640–10500 Å at resolving powers of 1100–1900. This broad wavelength coverage and integral-field capability make LRS2 well-suited for diverse astrophysical investigations, including studies of Lyman-alpha emitters (LAEs), planetary nebulae, stellar populations, and transient phenomena such as supernovae and optical counterparts to gravitational-wave events. LRS2 also supports follow-up work for large-scale surveys such as HETDEX \citep{Gebhardt2021} and enables research on high-redshift galaxies, active galactic nuclei (AGN), brown dwarfs, and local emission-line galaxies.

Panacea automates the complete end-to-end reduction of raw LRS2 CCD frames into science-ready spectra and data cubes. It executes daily on the Texas Advanced Computing Center (TACC) systems and supports reproducible user-driven reductions, providing the community with uniform and verified data products.

# Statement of Need

Modern multi-arm integral-field spectrographs such as LRS2 generate thousands of spectra per night across four independent channels, each requiring distinct calibration and extraction strategies. Manual reduction of such data is impractical for routine operations. **Panacea** provides a reliable, standardized, and fully automated framework that performs all calibration, extraction, and combination steps to produce consistent, reproducible results.

The primary users of Panacea are **LRS2 observers and other Hobby–Eberly Telescope researchers**, but astronomers working with other IFU spectrographs may also individual algorithms useful for adapting to similar instruments. Since its introduction on **1 January 2019**, Panacea has continuously processed all production LRS2 data on TACC and has supported more than **57 refereed publications** citing Panacea products from the Hobby–Eberly Telescope as of **October 2025**, underscoring its maturity and broad scientific impact.

# State of the Field

Several community pipelines address spectroscopic data reduction, including **PypeIt** \citep{Prochaska2020}, a flexible long-slit and echelle reduction framework, and **PyWiFeS**, designed for the WiFeS integral-field spectrograph. **Remedy** \citep{Zeimann2024} provides a related system for the VIRUS spectrographs on HET. Panacea differs from these packages by integrating both CCD-level calibration and fiber-based spectral reconstruction in a single framework, specifically tailored for the **dual-arm, multi-amplifier IFU design of LRS2**. Its automation and provenance tracking are optimized for nightly survey operations as well as principal-investigator data reduction.

# Methods and Implementation

Panacea is written in **Python 3** and orchestrates each stage of the LRS2 reduction sequence, from overscan and bias subtraction to flux calibration, within a single, reproducible workflow. It models and removes amplifier-dependent electronic offsets, traces and extracts fiber profiles using an optimal extraction algorithm, and derives wavelength solutions from arc-lamp exposures with sub-pixel precision. Flat-fielding and fiber-to-fiber normalization correct for spatial and throughput variations, while an empirical two-dimensional sky model minimizes residuals from bright sky lines. The resulting spectra are relatively flux-calibrated using standard response curves and put on absolute scale using guider-based transparency estimates as well as mirror illumination models of the fixed altitude HET. 

Panacea produces multi-extension FITS files that include extracted spectra, sky models, error frames, and ancillary diagnostic information. It automatically detects targets in the IFU frame across the channel wavelengths and extracts a spectrum using an optimal weighting algorithm if the source has a signal to noise greater than 5.  Finally, we also provide data cubes for each channel.  These data products enable both scientific analysis and verification of calibration quality.

# Validation

The accuracy of Panacea reductions is assessed through repeat observations of spectrophotometric standard stars and cross-channel consistency tests. Validation confirms stable wavelength and flux calibration performance across all four spectrograph arms, demonstrating that Panacea delivers reliable, science-quality data for both nightly and archival use.

# Author Contributions

- **Conceptualization, Software, Validation, Writing – Original Draft:** Greg Zeimann

# Acknowledgements

The author thanks the **Hobby–Eberly Telescope** operations staff, the **HETDEX collaboration**, and administrators at the **Texas Advanced Computing Center (TACC)** for their essential support of the automated reduction framework. Additional gratitude is extended to HET instrument scientists and LRS2 observers for providing continuous feedback that shaped Panacea’s development.

# References
See `paper.bib`.
