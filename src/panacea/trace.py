"""Tracing utilities.

Migrated functions per function_map.md.
"""
from __future__ import annotations

import glob
import os.path as op
from datetime import datetime

import numpy as np


def get_trace_shift(sci_array: np.ndarray, flat: np.ndarray, array_trace: np.ndarray, Yx: np.ndarray) -> np.ndarray:
    """Compute per-row shifts aligning science traces to flat-field traces.

    Uses a quadratic interpolation around the central pixel to estimate peak
    positions of both science and flat images around each trace row, then
    computes median offset and fits a quadratic relation vs. median trace
    position to evaluate shifts at coordinates Yx.

    Args:
        sci_array: 2D science image.
        flat: 2D flat-field image.
        array_trace: 2D trace positions (fibers x columns).
        Yx: 1D array of y-indices where shifts should be evaluated.

    Returns:
        1D array of shifts evaluated at Yx.
    """
    YM, XM = np.indices(flat.shape)
    inds = np.zeros((3, array_trace.shape[0], array_trace.shape[1]))
    XN = np.round(array_trace)
    inds[0] = XN - 1.0
    inds[1] = XN + 0.0
    inds[2] = XN + 1.0
    inds = np.array(inds, dtype=int)
    Trace = array_trace * 0.0
    FlatTrace = array_trace * 0.0
    N = YM.max()
    x = np.arange(array_trace.shape[1])
    for i in np.arange(Trace.shape[0]):
        sel = YM[inds[0, i, :], x] >= 0.0
        sel = sel * (YM[inds[2, i, :], x] < N)
        xmax = (
            YM[inds[1, i, sel], x[sel]]
            - (sci_array[inds[2, i, sel], x[sel]] - sci_array[inds[0, i, sel], x[sel]])
            / (2.0 * (sci_array[inds[2, i, sel], x[sel]] - 2.0 * sci_array[inds[1, i, sel], x[sel]] + sci_array[inds[0, i, sel], x[sel]]))
        )
        Trace[i, sel] = xmax
        xmax = (
            YM[inds[1, i, sel], x[sel]]
            - (flat[inds[2, i, sel], x[sel]] - flat[inds[0, i, sel], x[sel]])
            / (2.0 * (flat[inds[2, i, sel], x[sel]] - 2.0 * flat[inds[1, i, sel], x[sel]] + flat[inds[0, i, sel], x[sel]]))
        )
        FlatTrace[i, sel] = xmax
    mid = int(Trace.shape[1] / 2)
    shifts = np.nanmedian((FlatTrace - Trace)[:, mid - 200 : mid + 200], axis=1)
    shifts = np.polyval(np.polyfit(np.nanmedian(FlatTrace, axis=1), shifts, 2), Yx)
    return shifts



def get_trace_reference(specid: str, ifuslot: str, ifuid: str, amp: str, obsdate: str, virusconfig: str = '/work/03946/hetdex/maverick/virus_config'):
    """Locate and load the closest-in-time reference trace file.

    Args:
        specid: Spectrograph ID.
        ifuslot: IFU slot.
        ifuid: IFU identifier.
        amp: Amplifier identifier.
        obsdate: Observation date as YYYYMMDD.
        virusconfig: Base path to virus configuration.

    Returns:
        ndarray with reference trace rows: columns [col, flag] per fiber.
    """
    files = glob.glob(op.join(virusconfig, 'Fiber_Locations', '*', f'fiber_loc_{specid}_{ifuslot}_{ifuid}_{amp}.txt'))
    dates = [op.basename(op.dirname(fn)) for fn in files]
    obsdate_dt = datetime(int(obsdate[:4]), int(obsdate[4:6]), int(obsdate[6:]))
    timediff = np.zeros((len(dates),))
    for i, datei in enumerate(dates):
        d = datetime(int(datei[:4]), int(datei[4:6]), int(datei[6:]))
        timediff[i] = abs((obsdate_dt - d).days)
    ref_file = np.loadtxt(files[int(np.argmin(timediff))])
    return ref_file


def get_trace(twilight: np.ndarray, specid: str, ifuslot: str, ifuid: str, amp: str, obsdate: str):
    """Compute per-fiber trace positions across detector columns.

    Args:
        twilight: 2D twilight flat image.
        specid: Spectrograph ID.
        ifuslot: IFU slot.
        ifuid: IFU identifier.
        amp: Amplifier identifier.
        obsdate: Observation date as YYYYMMDD.

    Returns:
        Tuple of (trace, ref) where trace is (fibers x columns) and ref is the
        reference array used to seed missing fibers.
    """
    import numpy as np

    ref = get_trace_reference(specid, ifuslot, ifuid, amp, obsdate)
    N1 = int((ref[:, 1] == 0.0).sum())
    good = np.where(ref[:, 1] == 0.0)[0]

    def _get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.0
        inds[1] = XN + 0.0
        inds[2] = XN + 1.0
        inds = np.array(inds, dtype=int)
        Trace = (
            YM[inds[1]]
            - (flat[inds[2]] - flat[inds[0]])
            / (2.0 * (flat[inds[2]] - 2.0 * flat[inds[1]] + flat[inds[0]]))
        )
        return Trace

    image = twilight
    N = 40
    xchunks = np.array([np.mean(x) for x in np.array_split(np.arange(image.shape[1]), N)])
    chunks = np.array_split(image, N, axis=1)
    flats = [np.median(chunk, axis=1) for chunk in chunks]
    Trace = np.zeros((len(ref), len(chunks)))
    k = 0
    for flat, x in zip(flats, xchunks):
        diff_array = flat[1:] - flat[:-1]
        loc = np.where((diff_array[:-1] > 0.0) * (diff_array[1:] < 0.0))[0]
        peaks = flat[loc + 1]
        loc = loc[peaks > 0.1 * np.median(peaks)] + 1
        trace = _get_trace_chunk(flat, loc)
        T = np.zeros((len(ref)))
        if len(trace) == N1:
            T[good] = trace
            for missing in np.where(ref[:, 1] == 1)[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = T[good[gind]] + ref[missing, 0] - ref[good[gind], 0]
        Trace[:, k] = T
        k += 1
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] > 0.0
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)
    return trace, ref
