"""
Deterministic tiny synthetic LRS2-like frames for tests.

We generate a minimal set of FITS images for a single channel with realistic
headers and tiny dimensions to keep the repo small and tests fast.

Generated frames:
- bias: 128x64, constant level with small 2D structure
- flat: continuum illumination with fiber-profile stripes
- arc: sparse emission lines along dispersion
- twilight: smooth continuum + weak sky residual
- science: point source mapped across fiber profiles + sky background

The goal is to exercise I/O and light-weight parts of the pipeline without
shipping large binaries. These are not physically accurate, but they are
shape/metadata compatible for smoke tests.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from astropy.io import fits

RNG = np.random.default_rng(12345)

DEFAULT_SIZE = (64, 128)  # (ny, nx) keep small but non-trivial aspect


def _base_header(channel: str = "LRS2-R") -> fits.Header:
    hdr = fits.Header()
    hdr["INSTRUME"] = ("LRS2", "Instrument")
    hdr["CHANNEL"] = (channel, "LRS2 channel")
    hdr["DATE-OBS"] = ("2025-01-01T00:00:00", "Observation date")
    hdr["EXPTIME"] = (10.0, "Exposure time [s]")
    hdr["GAIN"] = (1.8, "e-/ADU")
    hdr["RDNOISE"] = (3.5, "e- rms")
    hdr["BSCALE"] = (1.0, "scale factor for data")
    hdr["BZERO"] = (0.0, "offset for data")
    # Minimal WCS-esque keywords for dispersion direction (x)
    hdr["CTYPE1"] = ("LINEAR", "Dispersion axis")
    hdr["CRPIX1"] = 1.0
    hdr["CRVAL1"] = 6500.0
    hdr["CDELT1"] = 1.5
    return hdr


def _fiber_traces(
    ny: int, nx: int, n_fibers: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate simple near-horizontal Gaussian fiber profiles.

    Returns (y_centers, sigma) arrays per fiber.
    """
    y_centers = np.linspace(5, ny - 6, n_fibers)
    sigma = np.full(n_fibers, 1.2)
    return y_centers, sigma


def _render_fibers(
    ny: int, nx: int, amp: float = 1000.0, n_fibers: int = 10
) -> np.ndarray:
    y = np.arange(ny)[:, None]
    y0, sig = _fiber_traces(ny, nx, n_fibers)
    img = np.zeros((ny, nx), dtype=np.float32)
    # Mild throughput variation per fiber
    thru = 0.8 + 0.4 * np.sin(np.linspace(0, 3.14, n_fibers))
    for i, (yc, s) in enumerate(zip(y0, sig)):
        prof = np.exp(-0.5 * ((y - yc) / s) ** 2).astype(np.float32)
        img += (amp * thru[i]) * prof
    return img


def _add_lines(
    img: np.ndarray, wavelengths: np.ndarray, hdr: fits.Header, strength: float = 2000.0
) -> None:
    nx = img.shape[1]
    crval = float(hdr.get("CRVAL1", 6500.0))
    cdelt = float(hdr.get("CDELT1", 1.5))
    # Map wavelengths to x indices
    x = np.round((wavelengths - crval) / cdelt).astype(int)
    for xi in x:
        if 0 <= xi < nx:
            img[:, xi] += strength


def _bias_frame(shape: Tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    ygrid, xgrid = np.mgrid[0:ny, 0:nx]
    plane = 500.0 + 0.3 * (xgrid / nx) + 0.2 * (ygrid / ny)
    noise = RNG.normal(0, 2.0, size=shape)
    return (plane + noise).astype(np.float32)


def _flat_frame(shape: Tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    img = 2000.0 * np.ones((ny, nx), dtype=np.float32)
    img += 100.0 * (np.linspace(0, 1, nx, dtype=np.float32)[None, :])
    img += _render_fibers(ny, nx, amp=500.0)
    img *= RNG.normal(1.0, 0.01, size=(ny, nx)).astype(np.float32)
    return img


def _arc_frame(shape: Tuple[int, int], hdr: fits.Header) -> np.ndarray:
    ny, nx = shape
    img = np.zeros((ny, nx), dtype=np.float32)
    img += _render_fibers(ny, nx, amp=50.0)
    # A few bright lines spanning the field
    lines = np.array([6200.0, 6500.0, 6800.0, 7000.0])
    _add_lines(img, lines, hdr, strength=1500.0)
    # small noise
    img += RNG.normal(0, 3.0, size=(ny, nx)).astype(np.float32)
    return img


def _twilight_frame(shape: Tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    base = 1000.0 * np.ones((ny, nx), dtype=np.float32)
    base += 30.0 * np.sin(np.linspace(0, 10, nx, dtype=np.float32))[None, :]
    base += _render_fibers(ny, nx, amp=150.0)
    base += RNG.normal(0, 2.0, size=(ny, nx)).astype(np.float32)
    return base


def _science_frame(shape: Tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    img = 200.0 * np.ones((ny, nx), dtype=np.float32)  # sky background
    img += _render_fibers(ny, nx, amp=120.0)
    # Point source across a few columns centered in x
    y = np.arange(ny)[:, None]
    x = np.arange(nx)[None, :]
    yc = ny / 2 - 2.0
    xc = nx / 2
    psf_y = np.exp(-0.5 * ((y - yc) / 1.8) ** 2)
    psf_x = np.exp(-0.5 * ((x - xc) / 3.0) ** 2)
    img += (800.0 * psf_y * psf_x).astype(np.float32)
    img += RNG.normal(0, 2.0, size=(ny, nx)).astype(np.float32)
    return img


def write_sample_dataset(
    outdir: Path, channel: str = "LRS2-R", shape: Tuple[int, int] = DEFAULT_SIZE
) -> Dict[str, Path]:
    """Create small synthetic frames into outdir. Returns dict of paths.

    Files created (FITS, float32 primary):
    - bias.fits
    - flat.fits
    - arc.fits
    - twilight.fits
    - science.fits
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hdr = _base_header(channel)

    paths: Dict[str, Path] = {}

    def _write(name: str, data: np.ndarray):
        h = hdr.copy()
        h["IMAGETYP"] = name
        fits.PrimaryHDU(data=data.astype(np.float32), header=h).writeto(
            outdir / f"{name}.fits", overwrite=True
        )
        paths[name] = outdir / f"{name}.fits"

    _write("bias", _bias_frame(shape))
    _write("flat", _flat_frame(shape))
    _write("arc", _arc_frame(shape, hdr))
    _write("twilight", _twilight_frame(shape))
    _write("science", _science_frame(shape))

    return paths


__all__ = ["write_sample_dataset", "DEFAULT_SIZE"]
