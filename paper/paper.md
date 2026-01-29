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
date: 2026-01-29
bibliography: paper.bib
---

# Summary

Panacea is the automated data-reduction pipeline for the Low Resolution Spectrograph 2 (LRS2; @Chonis:2016) on the Hobby–Eberly Telescope (HET; @Ramsey:1998; @Hill:2021). LRS2 is a fiber-fed integral-field spectrograph mounted on the 10-m HET, designed to deliver spatially resolved optical spectroscopy across a broad wavelength range.

The instrument’s four spectrograph channels (UV, Orange, Red, and Far Red) span 3640–10500 Å at resolving powers of approximately 1100–1900. This configuration supports a wide range of astrophysical applications, including studies of Lyman-alpha emitters, planetary nebulae, stellar populations, and time-domain phenomena such as supernovae and optical counterparts to gravitational-wave sources. LRS2 also provides follow-up spectroscopy for large surveys such as HETDEX (@Gebhardt:2021) and enables detailed investigations of high-redshift galaxies, active galactic nuclei, brown dwarfs, and emission-line systems in the nearby Universe.

Panacea automates the transformation of raw LRS2 CCD frames into science-ready spectra and data cubes. The pipeline executes daily on Texas Advanced Computing Center (TACC) systems, producing uniform, provenance-tracked data products that are immediately suitable for scientific analysis.

# Statement of Need

Modern multi-arm integral-field spectrographs generate thousands of spectra per night across multiple detectors, each requiring careful calibration, extraction, and combination. Manual or ad hoc reduction of such data is time-intensive and prone to inconsistencies, particularly in an operational observatory environment. Panacea addresses this challenge by providing a standardized, automated framework that performs all core reduction steps in a consistent and reproducible manner.

The primary users of Panacea are LRS2 observers and the broader Hobby–Eberly Telescope community. However, astronomers working with other fiber-fed IFU spectrographs may adapt Panacea’s modular algorithms for similar systems. Since its deployment in 2019, Panacea has processed all production LRS2 data at HET and has contributed to more than 50 refereed publications as of October 2025, demonstrating its robustness and sustained scientific impact.

# State of the Field

Data-reduction pipelines for integral-field spectroscopy range from instrument-specific observatory systems to flexible, community-driven frameworks.

- **MUSE (ESO) pipeline.** The European Southern Observatory’s official MUSE pipeline (@Weilbacher:2020) provides a comprehensive, end-to-end reduction system for the 24-IFU image-slicer instrument.

- **KCWI pipelines.** The Keck Cosmic Web Imager is supported by the KCWI Data Reduction Pipeline (@Morrissey:2018; @Neill:2023), a Python-based system distributed through the Keck Observatory Archive. It addresses slicer geometry and calibration, while source detection is typically handled in post-cube analysis software.

- **MaNGA (SDSS-IV) DRP.** The MaNGA pipeline (@Law:2016) processed thousands of fiber bundles per night, performing wavelength calibration, sky subtraction, flux calibration, and rectified cube assembly. Its architecture established a widely adopted model for large-survey IFU reductions.

- **SDSS-V Local Volume Mapper (LVM).** The LVM Data Analysis Pipeline (DAP; @Sanchez:2025) extends the MaNGA framework to parsec-scale mapping of the Milky Way and nearby galaxies, combining calibration, data fusion, and distributed analysis optimized for wide-field mosaics.

- **PypeIt.** PypeIt (@Prochaska:2020) is a flexible, general-purpose spectroscopic pipeline supporting long-slit, multi-slit, and echelle data. While its modular calibration models are broadly applicable, it is not optimized for fiber-fed IFUs or multi-amplifier architectures such as LRS2.

- **Remedy.** Remedy (@Zeimann:2024) is the production reduction system for the VIRUS spectrographs on HET. It is optimized for massively multiplexed fiber spectroscopy and survey-scale operation for the HET VIRUS Parallel Survey (HETVIPS; @Zeimann:2024). Panacea and Remedy share common design patterns for amplifier handling and data provenance, but differ in scope: Panacea targets dual-arm IFU observations and observer-driven data products, whereas Remedy is tailored for survey pipelines.

These pipelines span a range of use cases, from large survey operations to general-purpose spectroscopy. However, few systems are designed to support multi-arm, fiber-fed IFU instruments in an operational observatory setting with automated daily reductions and built-in source detection.

**Panacea’s distinction.**
1. **Integrated CCD-to-spectrum reduction for LRS2.** Panacea combines CCD-level calibration (bias subtraction, gain correction, tracing, wavelength calibration, and throughput normalization) with fiber-based extraction in a unified framework.
2. **Built-in automatic target detection on IFU data.** Unlike many pipelines that rely on post-cube detection, Panacea performs PSF- and fiber-aware automatic detection directly on IFU frames using matched filtering and optimal fiber weighting.
3. **Automated daily reductions with reproducibility.** Panacea executes automatically each morning on TACC systems, delivering consistent, validated spectra and data cubes with full provenance tracking.

# Software Design

Panacea is written in Python 3 and orchestrates each stage of the LRS2 reduction sequence independently for each spectrograph channel, from overscan and bias subtraction through flux calibration, within a single reproducible workflow. The pipeline models and removes amplifier-dependent offsets, traces and optimally extracts fiber profiles, and derives wavelength solutions from arc-lamp exposures with sub-pixel precision. Flat-fielding and fiber-to-fiber normalization correct throughput variations, while a two-dimensional sky model suppresses residuals from bright night-sky emission lines. The extracted spectra are relatively flux calibrated using standard response curves and placed on an absolute scale using guider-based transparency estimates and mirror illumination models appropriate for the fixed-altitude HET.

Panacea produces multi-extension FITS files containing extracted spectra, sky models, variance estimates, and diagnostic data products. For automatic target detection, the pipeline masks bright skylines and cosmic rays, smooths spectra spectrally, and constructs per-fiber signal-to-noise images to identify the most significant wavelength slice. It then collapses a narrow spectral window, fits a two-dimensional Gaussian to estimate source centroid and spatial extent, corrects for differential atmospheric refraction, and performs optimal extraction when the detection exceeds a signal-to-noise threshold of five. Rectified data cubes are also generated for each spectrograph channel to support science analysis and calibration verification.

# Quality Control

Pipeline performance is validated through repeated observations of spectrophotometric standard stars and cross-channel consistency tests. These assessments demonstrate stable wavelength and flux calibration across all four spectrograph arms, confirming that Panacea delivers reproducible, science-quality data products suitable for both nightly operations and archival analyses.

# Research Impact Statement

Panacea has been in continuous production use at the HET since 2019, reducing all LRS2 data as part of daily automated operations at the Texas Advanced Computing Center (TACC). This operational role has enabled immediate, uniform availability of calibrated spectra and cubes to observers and the HET community. As of October 2025, LRS2 data reduced with Panacea have contributed to more than 50 refereed publications, evidencing sustained scientific impact across diverse programs (e.g., transients, emission‑line galaxies, nearby nebulae, and survey follow‑up). The software is openly available under a permissive license with installation instructions, an audited dependency environment (environment.yml), command‑line entry points, API documentation, and a reproducible documentation build, signaling community readiness and enabling external validation and reuse.

# AI Usage Disclosure

Generative AI tools were used to assist with documentation, paper editing for this submission (e.g., revising prose), and with code re-factoring from the initial monolithic script built in 2019. All AI‑assisted text was reviewed, edited, and verified by the author for technical accuracy and correctness. No generative AI tools were used to write the scientific software itself; the Panacea codebase is human‑written and maintained, with behavior validated by routine operational checks, tests, and documentation builds.

# Acknowledgements

The author thanks the HET operations staff and the Texas Advanced Computing Center (TACC) for their support of the automated reduction infrastructure. Additional thanks go to HET instrument scientists and LRS2 observers whose feedback guided Panacea’s development and long-term refinement.

# References

See `paper.bib`. The archived software release corresponding to this paper is Panacea v1.0.1 (Zenodo concept DOI: 10.5281/zenodo.18250411).
