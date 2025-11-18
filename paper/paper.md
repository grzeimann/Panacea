---
title: 'Panacea: The LRS2 Data Reduction Pipeline for the Hobby–Eberly Telescope'
tags:
  - astronomy
  - spectroscopy
  - data reduction
  - IFU
  - Hobby–Eberly Telescope
authors:
  - name: Gregory Zeimann
    orcid: 0000-0003-2307-0629
    affiliation: 1
affiliations:
  - name: Hobby–Eberly Telescope, University of Texas, Austin
    index: 1
date: 2025-11-12
bibliography: paper.bib
---

# Summary

Panacea is the automated data-reduction pipeline for the Low Resolution Spectrograph 2 [LRS2; @Chonis:2016]  on the Hobby–Eberly Telescope (HET) [@Ramsey:1998; @Hill:2021]. LRS2 is a fiber-fed integral-field spectrograph mounted on the 10-m HET, designed to capture spatially resolved optical spectra across a broad wavelength range. 

The instrument’s four spectrograph channels (UV, Orange, Red, and Far Red) span 3640–10500 Å at resolving powers of 1100–1900. This configuration supports a range of astrophysical applications, including studies of Lyman-alpha emitters, planetary nebulae, stellar populations, and transient events such as supernovae and optical counterparts to gravitational-wave sources. LRS2 also provides follow-up spectroscopy for large-scale surveys such as HETDEX [@Gebhardt:2021] and enables studies of high-redshift galaxies, active galactic nuclei, brown dwarfs, and emission-line systems in the nearby Universe.

Panacea automates the reduction of raw LRS2 CCD frames into science-ready spectra and data cubes. It executes daily on the Texas Advanced Computing Center (TACC) systems providing uniform, verified data products for scientific analysis.

# Statement of Need

Modern multi-arm integral-field spectrographs generate thousands of spectra per night across multiple detectors, each requiring dedicated calibration and extraction. Manual reduction of such data is time-intensive and inconsistent for routine operations. Panacea provides a standardized, automated framework that performs calibration, extraction, and combination steps in a consistent, reproducible workflow.

The primary users of Panacea are LRS2 observers and other Hobby–Eberly Telescope researchers, but astronomers working with other IFU spectrographs may adapt its modular algorithms for similar systems. Since its deployment in 2019, Panacea has processed all production LRS2 data at HET and contributed to more than 50 refereed publications as of October 2025, demonstrating its reliability and broad use within the community.

# State of the Field

Modern IFU pipelines span observatory-maintained systems and general-purpose frameworks.

- **MUSE (ESO) pipeline.** The European Southern Observatory’s official MUSE pipeline [@Weilbacher:2020] provides an end-to-end reduction system for the 24-IFU image-slicer instrument. 

- **KCWI pipelines.** The Keck Cosmic Web Imager is supported by the KCWI Data Reduction Pipeline [@Morrissey:2018], implemented in Python and distributed via the Keck Observatory Archive, handles slicer geometry and calibration but generally defer source detection to post-cube analysis software.

- **MaNGA (SDSS-IV) DRP.** The MaNGA pipeline [@Law:2016] processed thousands of fiber bundles per night, performing wavelength calibration, sky subtraction, flux calibration, and rectified cube assembly. Its architecture established a model for large-survey IFU reductions.

- **SDSS-V Local Volume Mapper (LVM).** The LVM DRP and DAP [@Sanchez:2024; @Sanchez:2025] extend the MaNGA framework to parsec-scale mapping of the Milky Way and nearby galaxies, combining calibration, data fusion, and distributed analysis optimized for wide-field mosaics.

- **PypeIt.** PypeIt [@Prochaska:2020] is a flexible, general-purpose spectroscopic pipeline supporting long-slit, multi-slit, and echelle data. Its modular design and calibration models have broad applicability, though it is not optimized for fiber-fed IFUs or multi-amplifier architectures such as LRS2.

- **Remedy.** Remedy [@Zeimann:2024] is the production reduction system for the VIRUS spectrographs on HET. It is optimized for massively multiplexed fiber spectroscopy and survey-scale operation for the HET VIRUS Parallel Survey [HETVIPS; @Zeimann:2024], emphasizing throughput, automated calibration, and efficient sky modeling. Panacea and Remedy share common software patterns for amplifier handling and data provenance but differ in scope: Panacea focuses on dual-arm IFU observations and flexible, observer-driven data products.

**Panacea’s distinction.**
1. **Integrated CCD–to–spectrum reduction for LRS2.** Panacea combines CCD-level calibration (bias removal, gain and trace solutions, wavelength calibration, and throughput normalization) with fiber-based extraction in a unified framework.
2. **Built-in automatic target detection on IFU frames.** Unlike most pipelines that rely on post-cube detection, Panacea includes PSF- and fiber–aware automatic detection directly on the IFU frame, identifying and extracting sources using matched filtering and optimal fiber weighting.
3. **Automated daily reductions with reproducibility.** Panacea executes automatically each morning on TACC systems, producing consistent, provenance-tracked spectra and cubes that are ready for principal investigator use.

# Methods and Implementation

Panacea is written in Python 3 and orchestrates each stage of the LRS2 reduction sequence for each channel independently, from overscan and bias subtraction to flux calibration, within a single reproducible workflow. It models and removes amplifier-dependent offsets, traces and extracts fiber profiles with an optimal extraction algorithm, and derives wavelength solutions from arc-lamp exposures with sub-pixel precision. Flat-fielding and fiber-to-fiber normalization correct for throughput variations, while a two-dimensional sky model minimizes residuals from bright sky lines. The resulting spectra are relatively flux-calibrated using standard response curves and placed on an absolute scale using guider-based transparency estimates and mirror illumination models of the fixed-altitude HET. 

Panacea produces multi-extension FITS files that include extracted spectra, sky models, error frames, and diagnostic extensions. For automatic target detection in the IFU, the pipeline masks bright skylines and cosmic rays, smooths spectra with a Gaussian kernel, and constructs a per-fiber signal-to-noise image to locate the most significant wavelength slice away from edges. It collapses a narrow spectral window, fits a two-dimensional Gaussian to nearby fibers to estimate centroid and apparent size, recenters for differential atmospheric refraction, and performs an optimal extraction when the detection exceeds S/N > 5. Data cubes are also generated for each spectrograph channel, supporting both science analysis and calibration verification.

# Validation

Performance validation is conducted through repeat observations of spectrophotometric standard stars and cross-channel consistency tests. Results confirm stable wavelength and flux calibration across all four spectrograph arms, demonstrating that Panacea delivers reproducible, science-quality data for nightly and archival use.

# Acknowledgements

The author thanks the Hobby–Eberly Telescope operations staff and administrators at the Texas Advanced Computing Center (TACC) for their support of the automated reduction system. Additional thanks go to HET instrument scientists and LRS2 observers for their feedback that guided Panacea’s development.

# References
See `paper.bib`.
