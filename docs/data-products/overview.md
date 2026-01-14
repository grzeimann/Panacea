# Data Products

This page summarizes the main outputs produced by Panacea and their contents.

## Primary products
- spectrum*.fits — produced for all exposures and channels.
```text
row1: wavelength (air)
row2: extracted object spectrum (f_lambda: ergs/s/cm^2/A)
row3: extracted sky spectrum from same aperture and weighting as object (s_lambda: ergs/s/cm^2/A)
row4: error for extracted object spectrum (e_f_lambda: ergs/s/cm^2/A)
row5: error for extracted sky spectrum (e_s_lambda: ergs/s/cm^2/A)
row6: response function (ergs / e-)
```

- multi*{uv,orange,red,farred}.fits — multi-extension FITS containing:

```
Rectified Spectra: flux calibrated spectrum (object + sky) for each fiber
Rectified Sky Model:flux calibrated sky spectrum for each fiber
Rectified Sky Subtracted Spectra: flux calibrated sky subtracted spectrum for each fiber
Rectified Error Frame: flux calibrated error spectrum for each fiber
Collapsed image: a collapsed frame for visualization of the source(s)
Positions (IFU, Focal, Sky): ifu x and y positions, focal x and y position, and ra and dec
Extracted Spectra and Response: This is identical to the spectrum*.fits extension above
ADR: The atmospheric differential refraction as a function of wavelength.  The columns are wavelength, x_adr, y_adr
CCD Wavelength: The wavelength of each pixel in the 2d frame
Image: the initial reduction of the 2d raw frame.
Flat Fielded image: same as the image frame above but divided by the flat field (fiber profile and fiber to fiber normalization)
Central Trace Pixels: location of the pixels for each fiber (central two pixels)
Cosmics: identified cosmics in the central four pixels of the trace
Unrectified Spectra: Unrectified, uncalibrated spectra for each fiber
```

- *cube*.fits — per-channel spectral datacubes.
  - Spatial axes: IFU plane coordinates covering roughly 7×11 arcsec per LRS2 IFU (one cube per channel). The pixel scale is set by the reconstruction grid.
  - Spectral axis: wavelength in vacuum Angstroms (Å) as the third dimension.
  - Best use: quickly identify your target and strong emission features for quality control and visualization. These cubes are meant for inspection rather than precise spectrophotometry or cross‑channel combination.

See also: [TACC Overview](../tacc/overview.md), [Algorithms](../algorithms/overview.md)
