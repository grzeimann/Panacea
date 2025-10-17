"""General routines orchestrating reduction steps.

Functions migrated per function_map.md.
"""
from __future__ import annotations

import os
import os.path as op
import glob
import tarfile
from datetime import datetime

import numpy as np
from astropy.io import fits

from .utils import truncate_list


def get_ifucenfile(side: str, amp: str, virusconfig: str = "/work/03946/hetdex/maverick/virus_config", skiprows: int = 4) -> np.ndarray:
    """Load IFU center positions for a given side and amplifier.

    Args:
        side: One of {"uv", "orange", "red", "farred"}.
        amp: Amplifier identifier ("LL", "LU", "RL", "RU").
        virusconfig: Base path to the VIRUS configuration directory.
        skiprows: Number of header rows to skip in the mapping text file.

    Returns:
        Array of (x, y) fiber positions for the requested amplifier half.
    """
    file_dict = {
        "uv": "LRS2_B_UV_mapping.txt",
        "orange": "LRS2_B_OR_mapping.txt",
        "red": "LRS2_R_NR_mapping.txt",
        "farred": "LRS2_R_FR_mapping.txt",
    }

    ifucen = np.loadtxt(op.join(virusconfig, "IFUcen_files", file_dict[side]), usecols=[0, 1, 2], skiprows=skiprows)

    if amp == "LL":
        return ifucen[140:, 1:3][::-1, :]
    if amp == "LU":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RL":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RU":
        return ifucen[140:, 1:3][::-1, :]
    raise ValueError(f"Unknown amp: {amp}")



def get_mirror_illumination(fn=None, default=51.4e4):
    """Compute mirror illumination area using hetillum (external), fallback to default.

    Args:
        fn: FITS filename or file-like to read header geometry parameters.
        default: Default area (cm^2) to return on failure.

    Returns:
        Mirror illumination area in cm^2.
    """
    try:
        F = fits.open(fn)
        names = ['RHO_STRT', 'THE_STRT', 'PHI_STRT', 'X_STRT', 'Y_STRT']
        r, t, p, x, y = [F[0].header[name] for name in names]
        cmd = (
            '/work/03730/gregz/maverick/illum_lib/hetillum -p'
            f' -x "[{x:.4f},{y:.4f},{p:.4f}]" "[{0.042:.4f},{0.014:.4f}]" 256'
        )
        mirror_illum = float(os.popen(cmd).read().split('\n')[0])
        area = mirror_illum * default
    except Exception:
        area = default
    return area


def get_mirror_illumination_guider(fn, exptime, default=51.4e4, path='/work/03946/hetdex/maverick'):
    """Estimate mirror illumination using guider tarballs around exposure time.

    Args:
        fn: Science filename to infer datetime.
        exptime: Exposure time (s).
        default: Default area if guider data not found.
        path: Base path containing guider tarballs by date.
    """
    try:
        M = []
        # infer date directory
        f = op.basename(fn)
        DT = f.split('_')[0]
        y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
        d0 = datetime(y, m, d, h, mi, s)
        tarfolder = op.join(path, 'gc1', '*.tar')
        tarfolder = glob.glob(tarfolder)
        if len(tarfolder) == 0:
            return default
        T = tarfile.open(tarfolder[0], 'r')
        init_list = sorted([name for name in T.getnames() if name.endswith('.fits')])
        final_list = []
        for t in init_list:
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.0) * (p < exptime + 10.0):
                final_list.append(t)
        # truncate list to 5
        final_list = truncate_list(final_list)
        for f2 in final_list:
            fobj = T.extractfile(T.getmember(f2))
            M.append(get_mirror_illumination(fobj))
        M = np.array(M)
        sel = M != default
        if sel.sum() > 0.0:
            area = float(np.mean(M[sel]))
        else:
            area = default
        return area
    except Exception:
        return default


def get_throughput(fn, exptime, path='/work/03946/hetdex/maverick'):
    """Compute average guider TRANSPAR during exposure window (two guiders).

    Returns a float, defaults to 1.0 if no active guider frames with TRANSPAR.
    """
    attr = ['GUIDLOOP', 'MJD', 'TRANSPAR']
    M = []
    f = op.basename(fn)
    DT = f.split('_')[0]
    y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
    d0 = datetime(y, m, d, h, mi, s)
    for gp in ['gc1', 'gc2']:
        tarfolder = op.join(path, gp, f'{gp}.tar')
        T = tarfile.open(tarfolder, 'r')
        init_list = sorted([name for name in T.getnames() if name.endswith('.fits')])
        final_list = []
        for t in init_list:
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.0) * (p < exptime + 10.0):
                final_list.append(t)
        final_list = truncate_list(final_list)
        for fnt in final_list:
            fobj = T.extractfile(T.getmember(fnt))
            f1 = fits.open(fobj)
            if f1[1].header.get('GUIDLOOP') == 'ACTIVE':
                M.append([])
                for att in attr:
                    M[-1].append(f1[1].header.get(att, 0))
    throughput = np.zeros((len(M),)) if len(M) else np.zeros((1,))
    for i, mi in enumerate(M):
        if len(mi) < 2:
            continue
        if mi[2] > 0.0:
            throughput[i] = mi[2]
    t = float(np.mean(throughput[throughput > 0.0])) if throughput.size else 1.0
    if not np.isfinite(t) or t > 1.1 or t < 0.1:
        t = 1.0
    return t
