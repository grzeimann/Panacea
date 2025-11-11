"""I/O utilities for Panacea.

Functions gradually migrated from full_lrs2_reduction.py per function_map.md.
"""
from __future__ import annotations

import os
import os.path as op
from datetime import datetime as dt, timedelta as td
import numpy as np
import glob as glob
import tarfile as tarfile
import fnmatch as fnmatch
import requests as requests
from astropy.io.votable import parse_single_table
from astropy.io import fits
import uuid


def get_tarname_from_filename(filename):
    """Return the associated tar archive path for a raw FITS path.

    The heuristic matches the legacy layout where raw files are stored in a
    tarball three directories above the file. For example, for a path like::

        /a/b/c/d/e/f.fits  ->  /a/b/c.tar

    Args:
        filename: Path to a FITS file within the raw tree.

    Returns:
        The derived .tar path.
    """
    tarname = op.dirname(op.dirname(op.dirname(filename))) + ".tar"
    return tarname



def get_filenames_from_tarfolder(tarfolder, path):
    """List member paths in the tar that match the given pattern.

    Args:
        tarfolder: Path to the .tar archive.
        path: Example path whose basename pattern will be matched inside the tar.

    Returns:
        Sorted list of pseudo-paths pointing to matched members.
    """


    T = tarfile.open(tarfolder, "r")
    names = T.getnames()
    matches = fnmatch.filter(names, op.join("*", op.basename(path)))
    matches = [op.join(op.dirname(tarfolder), match) for match in matches]
    matches = sorted(matches)
    T.close()
    return matches


def get_cal_path(pathname, date, ndays = 31):
    """Find calibration files around a date window, respecting controller swap.

    Args:
        pathname: Glob-like path template containing the date string.
        date: Base date as YYYYMMDD.
        ndays: Initial search window size.

    Returns:
        Sorted list of matching file paths extracted from tar archives.
    """


    date_ = dt(int(date[:4]), int(date[4:6]), int(date[6:]))
    date_controller_swap = dt(2024, 7, 22)
    flag_new = date_ > date_controller_swap
    filenames = []
    while len(filenames) == 0:
        datel = date_ - td(days=int(ndays / 2))
        for i in np.arange(ndays):
            ndate = datel + td(days=int(i))
            if flag_new and (ndate <= date_controller_swap):
                continue
            daten = f"{ndate.year:04d}{ndate.month:02d}{ndate.day:02d}"
            npath = pathname.replace(date, daten)
            tarpath = get_tarname_from_filename(npath)
            for tarname in sorted(glob.glob(tarpath)):
                filenames.append(get_filenames_from_tarfolder(tarname, npath))
        flat_list = [item for sublist in filenames for item in sublist]
        filenames = sorted(flat_list)
        ndays += 1
    return filenames


def create_image_header(wave, xgrid, ygrid, zgrid, func=None):
    """Create an ImageHDU-like object with pixel WCS for a 2D image.

    Args:
        wave: Unused (kept for compatibility with legacy signature).
        xgrid: 2D array of x pixel coordinates.
        ygrid: 2D array of y pixel coordinates.
        zgrid: 2D image to store in the HDU.
        func: HDU class to instantiate (default astropy.io.fits.ImageHDU).

    Returns:
        ImageHDU with CRVAL/CRPIX/CDELT keywords set.
    """
    if func is None:
        func = fits.ImageHDU
    hdu = func(np.array(zgrid, dtype="float32"))
    hdu.header["CRVAL1"] = xgrid[0, 0]
    hdu.header["CRVAL2"] = ygrid[0, 0]
    hdu.header["CRPIX1"] = 1
    hdu.header["CRPIX2"] = 1
    hdu.header["CTYPE1"] = "pixel"
    hdu.header["CTYPE2"] = "pixel"
    hdu.header["CDELT1"] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header["CDELT2"] = ygrid[1, 0] - ygrid[0, 0]
    return hdu


def write_cube(wave, xgrid, ygrid, zgrid, outname, he):
    """Write a 3D spectral cube FITS with basic header keywords.

    Args:
        wave: 1D wavelength array.
        xgrid: 2D x pixel grid.
        ygrid: 2D y pixel grid.
        zgrid: 3D data cube (nw, ny, nx) or equivalent.
        outname: Output FITS filename.
        he: Reference header (dict-like) to copy extra keywords from.
    """
    hdu = fits.PrimaryHDU(np.array(zgrid, dtype="float32"))
    hdu.header["CRVAL1"] = xgrid[0, 0]
    hdu.header["CRVAL2"] = ygrid[0, 0]
    hdu.header["CRVAL3"] = wave[0]
    hdu.header["CRPIX1"] = 1
    hdu.header["CRPIX2"] = 1
    hdu.header["CRPIX3"] = 1
    hdu.header["CTYPE1"] = "pixel"
    hdu.header["CTYPE2"] = "pixel"
    hdu.header["CTYPE3"] = "pixel"
    hdu.header["CDELT1"] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header["CDELT2"] = ygrid[1, 0] - ygrid[0, 0]
    hdu.header["CDELT3"] = wave[1] - wave[0]
    for key in he.keys():
        if key in hdu.header:
            continue
        if ("CCDSEC" in key) or ("DATASEC" in key):
            continue
        if ("BSCALE" in key) or ("BZERO" in key):
            continue
        try:
            hdu.header[key] = he[key]
        except Exception:
            continue
    hdu.writeto(outname, overwrite=True)


def panstarrs_query(ra_deg, dec_deg, rad_deg, mindet=1, maxsources=30000,
                    server=("https://archive.stsci.edu/panstarrs/search.php")):
    """Query Pan-STARRS DR1 @ MAST and return an Astropy Table.

    Args:
        ra_deg: Right ascension in degrees.
        dec_deg: Declination in degrees.
        rad_deg: Search radius in degrees.
        mindet: Minimum number of detections.
        maxsources: Maximum returned sources.
        server: Query URL.

    Returns:
        astropy.table.Table with results.
    """


    r = requests.get(
        server,
        params={
            "RA": ra_deg,
            "DEC": dec_deg,
            "SR": rad_deg,
            "max_records": maxsources,
            "outputformat": "VOTable",
            "ndetections": (">%d" % mindet),
        },
    )

    name = str(uuid.uuid4()) + ".xml"
    with open(name, "w") as outf:
        outf.write(r.text)

    data = parse_single_table(name)
    os.remove(name)
    return data.to_table(use_names_over_ids=True)
