"""Fit and plot the LRS2 UV and Orange dichroic correction."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar


def _load_lrs2_spectrum(obj):
    """
    Load a Panacea LRS2 spectrum product.

    Parameters
    ----------
    obj : str or object
        Either a filename or an object with a ``filename`` attribute.

    Returns
    -------
    dict
        Dictionary containing filename, wavelength, science spectrum,
        sky spectrum, and response function.

    Notes
    -----
    The assumed spectrum layout is:

    ``data[0]``
        Wavelength.
    ``data[1]``
        Calibrated science spectrum.
    ``data[2]``
        Calibrated sky spectrum.
    ``data[-1]``
        Response function.
    """
    if hasattr(obj, "filename"):
        filename = obj.filename
    else:
        filename = str(obj)

    filename = filename.replace("multi", "spectrum")

    with fits.open(filename) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)

    wave = data[0].copy()
    science = data[1].copy()
    sky = data[2].copy()
    response = data[-1].copy()

    order = np.argsort(wave)

    return {
        "filename": filename,
        "wave": wave[order],
        "science": science[order],
        "sky": sky[order],
        "response": response[order],
    }


def _make_interp(wave, values):
    """
    Build a shape-preserving interpolator from finite samples.

    Parameters
    ----------
    wave : array_like
        Wavelength coordinates.
    values : array_like
        Values sampled at ``wave``.

    Returns
    -------
    scipy.interpolate.PchipInterpolator
        Shape-preserving cubic interpolator with extrapolation.
    """
    good = np.isfinite(wave) & np.isfinite(values)

    x = wave[good]
    y = values[good]

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    keep = np.concatenate(([True], np.diff(x) > 0))
    x = x[keep]
    y = y[keep]

    return PchipInterpolator(x, y, extrapolate=True)


def _repair_orange_response(wave, response, enabled=False):
    """Optionally repair the blue Orange response edge with a quadratic fit."""
    response_original = np.asarray(response, dtype=float).copy()
    response_repaired = response_original.copy()

    if not enabled:
        return response_original, response_repaired

    trusted = (wave >= 4645.0) & (wave <= 4655.0) & np.isfinite(response_original)
    if np.count_nonzero(trusted) < 3:
        raise ValueError(
            "At least three finite Orange response samples are required from 4645--4655 Å."
        )

    centered_wave = wave - 4650.0
    coefficients = np.polyfit(centered_wave[trusted], response_original[trusted], 2)
    repair_mask = np.isfinite(wave) & (wave < 4645.0)
    response_repaired[repair_mask] = np.polyval(
        coefficients, centered_wave[repair_mask]
    )

    return response_original, response_repaired


def fit_lrs2_dichroic(
    uv,
    orange,
    fit_range=(4635.0, 4645.0),
    orange_fit_min=4638.0,
    shift_bounds=(-3.0, 3.0),
    scale_bounds=(0.95, 1.05),
    throughput_floor=0.02,
    min_valid_fraction=0.6,
    apply_grey_scale=True,
    scan_points=121,
    repair_orange_response=False,
):
    """
    Fit a wavelength displacement of the LRS2-B dichroic response.

    The UV and orange response curves are shifted by the same wavelength
    displacement. The fit is determined from agreement between the sky
    spectra in the channel overlap.

    The response model is

    ``R_new(lambda) = R_nominal(lambda - delta_lambda)``.

    If the existing spectrum has already been calibrated using
    ``R_nominal``, the corrected spectrum is

    ``F_new = F_old * R_new / R_nominal``.

    At every trial value of ``delta_lambda``, the relative orange sky
    normalization is solved analytically. After the sky-derived shift is
    fixed, a separate science normalization is solved analytically.

    Parameters
    ----------
    uv, orange : str or object
        Spectrum filenames, multi filenames, or objects with a
        ``filename`` attribute.
    fit_range : tuple of float, optional
        Wavelength interval used to compare the channels, in Angstrom.
    orange_fit_min : float, optional
        Lower wavelength cutoff for Orange samples used in the fit, in
        Angstrom. This defaults to 4638 Angstrom because the Orange
        spectrum is treated as unreliable below that wavelength.
    shift_bounds : tuple of float, optional
        Allowed dichroic wavelength displacement, in Angstrom.
    scale_bounds : tuple of float, optional
        Allowed multiplicative orange-to-UV normalization.
    throughput_floor : float, optional
        Minimum relative throughput allowed in the fit.
    min_valid_fraction : float, optional
        Minimum fraction of the nominal overlap that must remain usable
        at a trial wavelength shift.
    apply_grey_scale : bool, optional
        Apply the fitted sky and science orange normalizations to the
        corresponding returned orange spectra and science response.
    scan_points : int, optional
        Number of shifts evaluated for the diagnostic cost curve.
    repair_orange_response : bool, optional
        Fit a quadratic to the original Orange response from 4645 to
        4655 Angstrom and use its blueward extrapolation below 4645
        Angstrom during the fit. The default leaves the response unchanged.

    Returns
    -------
    dict
        Dictionary containing corrected UV and orange products, fit
        parameters, fit diagnostics, and the cost-function scan.
    """
    uv_data = _load_lrs2_spectrum(uv)
    orange_data = _load_lrs2_spectrum(orange)

    wave_uv = uv_data["wave"]
    wave_orange = orange_data["wave"]

    response_uv_native = uv_data["response"]
    response_orange_original, response_orange_repaired = _repair_orange_response(
        wave_orange, orange_data["response"], enabled=repair_orange_response
    )

    response_uv = _make_interp(wave_uv, response_uv_native)
    response_orange_original_interp = _make_interp(
        wave_orange, response_orange_original
    )
    response_orange = _make_interp(wave_orange, response_orange_repaired)

    sky_uv = _make_interp(wave_uv, uv_data["sky"])
    sky_orange = _make_interp(wave_orange, orange_data["sky"])
    science_uv = _make_interp(wave_uv, uv_data["science"])
    science_orange = _make_interp(wave_orange, orange_data["science"])

    overlap_min = max(fit_range[0], np.nanmin(wave_uv), np.nanmin(wave_orange))
    overlap_max = min(fit_range[1], np.nanmax(wave_uv), np.nanmax(wave_orange))

    if overlap_max <= overlap_min:
        raise ValueError(
            "No nominal wavelength overlap. "
            f"Requested range is {fit_range}; "
            f"available range is {overlap_min:.2f}--"
            f"{overlap_max:.2f}."
        )

    step_uv = np.nanmedian(np.abs(np.diff(wave_uv)))
    step_orange = np.nanmedian(np.abs(np.diff(wave_orange)))
    wavelength_step = max(step_uv, step_orange)

    grid = np.arange(overlap_min, overlap_max + 0.5 * wavelength_step, wavelength_step)

    sky_uv_nominal = sky_uv(grid)
    sky_orange_nominal = sky_orange(grid)

    response_uv_nominal = response_uv(grid)
    response_orange_original_nominal = response_orange_original_interp(grid)
    response_orange_nominal = response_orange(grid)

    throughput_uv = np.full_like(response_uv_nominal, np.nan)
    throughput_orange = np.full_like(response_orange_nominal, np.nan)

    good_uv = np.isfinite(response_uv_nominal) & (response_uv_nominal > 0)
    good_orange = np.isfinite(response_orange_nominal) & (response_orange_nominal > 0)

    throughput_uv[good_uv] = 1.0 / response_uv_nominal[good_uv]
    throughput_orange[good_orange] = 1.0 / response_orange_nominal[good_orange]

    if np.any(np.isfinite(throughput_uv)):
        throughput_uv /= np.nanmax(throughput_uv)

    if np.any(np.isfinite(throughput_orange)):
        throughput_orange /= np.nanmax(throughput_orange)

    base_mask = (
        np.isfinite(sky_uv_nominal)
        & np.isfinite(sky_orange_nominal)
        & np.isfinite(response_uv_nominal)
        & np.isfinite(response_orange_nominal)
        & (response_uv_nominal > 0)
        & (response_orange_nominal > 0)
        & np.isfinite(throughput_uv)
        & np.isfinite(throughput_orange)
        & (throughput_uv >= throughput_floor)
        & (throughput_orange >= throughput_floor)
        & (grid >= orange_fit_min)
    )

    n_base = np.sum(base_mask)

    if n_base < 5:
        raise RuntimeError(
            f"Only {n_base} usable overlap samples remain. "
            "Try lowering throughput_floor."
        )

    relative_response_uv = response_uv_nominal / np.nanmedian(
        response_uv_nominal[base_mask]
    )
    relative_response_orange = response_orange_nominal / np.nanmedian(
        response_orange_nominal[base_mask]
    )

    weight = 1.0 / (relative_response_uv**2 + relative_response_orange**2)
    weight[~np.isfinite(weight)] = 0.0

    flux_scale = np.nanmedian(
        np.abs(
            np.concatenate((sky_uv_nominal[base_mask], sky_orange_nominal[base_mask]))
        )
    )

    if not np.isfinite(flux_scale) or flux_scale <= 0:
        flux_scale = 1.0

    def evaluate_delta(delta):
        """Evaluate the sky agreement for one dichroic shift."""
        shifted_response_uv = response_uv(grid - delta)
        shifted_response_orange = response_orange(grid - delta)

        correction_uv = shifted_response_uv / response_uv_nominal
        correction_orange = shifted_response_orange / response_orange_nominal

        corrected_sky_uv = sky_uv_nominal * correction_uv
        corrected_sky_orange = sky_orange_nominal * correction_orange

        valid = (
            base_mask
            & np.isfinite(shifted_response_uv)
            & np.isfinite(shifted_response_orange)
            & (shifted_response_uv > 0)
            & (shifted_response_orange > 0)
            & np.isfinite(correction_uv)
            & np.isfinite(correction_orange)
            & np.isfinite(corrected_sky_uv)
            & np.isfinite(corrected_sky_orange)
            & (weight > 0)
        )

        n_valid = np.sum(valid)

        empty_result = {
            "cost": np.inf,
            "scale": np.nan,
            "valid": valid,
            "uv": corrected_sky_uv,
            "orange": corrected_sky_orange,
            "uv_correction": correction_uv,
            "orange_correction": correction_orange,
        }

        if n_valid < 3:
            return empty_result

        valid_fraction = n_valid / n_base

        if valid_fraction < min_valid_fraction:
            return empty_result

        local_weight = weight[valid]
        orange_values = corrected_sky_orange[valid]
        uv_values = corrected_sky_uv[valid]

        denominator = np.sum(local_weight * orange_values**2)

        if not np.isfinite(denominator) or denominator <= 0:
            return empty_result

        scale = np.sum(local_weight * orange_values * uv_values) / denominator

        scale = np.clip(scale, scale_bounds[0], scale_bounds[1])

        residual = uv_values - scale * orange_values

        cost = np.sum(local_weight * residual**2) / np.sum(local_weight) / flux_scale**2

        # Penalize solutions that discard part of the overlap.
        cost *= n_base / n_valid

        return {
            "cost": cost,
            "scale": scale,
            "valid": valid,
            "uv": corrected_sky_uv,
            "orange": corrected_sky_orange,
            "uv_correction": correction_uv,
            "orange_correction": correction_orange,
        }

    def objective(delta):
        """Scalar objective used by the bounded optimizer."""
        return evaluate_delta(delta)["cost"]

    optimization = minimize_scalar(
        objective, bounds=shift_bounds, method="bounded", options={"xatol": 1e-4}
    )

    delta_lambda = float(optimization.x)
    best_fit = evaluate_delta(delta_lambda)

    if not np.isfinite(best_fit["cost"]):
        raise RuntimeError("Dichroic fit failed to find a valid solution.")

    orange_scale = float(best_fit["scale"])

    no_shift_fit = evaluate_delta(0.0)
    orange_scale_no_shift = float(no_shift_fit["scale"])

    science_uv_nominal = science_uv(grid)
    science_orange_nominal = science_orange(grid)
    science_uv_values = science_uv_nominal * best_fit["uv_correction"]
    shifted_response_orange = response_orange(grid - delta_lambda)
    science_orange_correction = np.full_like(shifted_response_orange, np.nan)
    np.divide(
        shifted_response_orange,
        response_orange_original_nominal,
        out=science_orange_correction,
        where=response_orange_original_nominal != 0,
    )
    science_orange_values = science_orange_nominal * science_orange_correction
    science_valid = (
        best_fit["valid"]
        & np.isfinite(science_uv_values)
        & np.isfinite(science_orange_values)
        & np.isfinite(science_orange_correction)
        & (science_orange_correction > 0)
    )

    if np.sum(science_valid) < 3:
        raise RuntimeError(
            "Science grey-scale fit failed to find enough valid overlap samples."
        )

    science_weight = weight[science_valid]
    science_orange_values = science_orange_values[science_valid]
    science_uv_values = science_uv_values[science_valid]
    science_denominator = np.sum(science_weight * science_orange_values**2)

    if not np.isfinite(science_denominator) or science_denominator <= 0:
        raise RuntimeError("Science grey-scale fit failed to find a valid solution.")

    science_scale = (
        np.sum(science_weight * science_orange_values * science_uv_values)
        / science_denominator
    )
    science_scale = np.clip(science_scale, scale_bounds[0], scale_bounds[1])

    def correct_native(
        channel, response_interp, science_scale=1.0, sky_scale=1.0, response_scale=1.0
    ):
        """
        Apply the dichroic correction and separate scales on the native grid.
        """
        wave = channel["wave"]
        nominal_response = channel["response"]

        shifted_response = response_interp(wave - delta_lambda)

        valid = (
            np.isfinite(nominal_response)
            & (nominal_response > 0)
            & np.isfinite(shifted_response)
            & (shifted_response > 0)
        )

        correction = np.ones_like(nominal_response)

        correction[valid] = shifted_response[valid] / nominal_response[valid]

        science_correction = correction * science_scale
        sky_correction = correction * sky_scale
        corrected_science = channel["science"] * science_correction
        corrected_sky = channel["sky"] * sky_correction

        corrected_response = nominal_response.copy() * response_scale

        corrected_response[valid] = shifted_response[valid] * response_scale

        return {
            "filename": channel["filename"],
            "wave": wave.copy(),
            "science": corrected_science,
            "sky": corrected_sky,
            "response": corrected_response,
            "correction": science_correction,
            "science_correction": science_correction,
            "sky_correction": sky_correction,
            "response_shifted": shifted_response,
            "valid_correction": valid,
            "science_original": channel["science"].copy(),
            "sky_original": channel["sky"].copy(),
            "response_original": nominal_response.copy(),
        }

    uv_output = correct_native(uv_data, response_uv)

    if apply_grey_scale:
        orange_sky_scale = orange_scale
        orange_science_scale = science_scale
    else:
        orange_sky_scale = 1.0
        orange_science_scale = 1.0

    # The calibrated Orange spectrum used the original response. Keep it as
    # the denominator while evaluating the shifted repaired response above.
    orange_output = correct_native(
        orange_data,
        response_orange,
        science_scale=orange_science_scale,
        sky_scale=orange_sky_scale,
        response_scale=orange_science_scale,
    )
    orange_output["response_repaired"] = response_orange_repaired.copy()

    mask_before = no_shift_fit["valid"]
    mask_after = best_fit["valid"]

    uv_before = no_shift_fit["uv"]
    orange_before = orange_scale_no_shift * no_shift_fit["orange"]

    uv_after = best_fit["uv"]
    orange_after = orange_scale * best_fit["orange"]

    residual_before = uv_before - orange_before
    residual_after = uv_after - orange_after

    rms_before = np.sqrt(np.nanmean(residual_before[mask_before] ** 2))
    rms_after = np.sqrt(np.nanmean(residual_after[mask_after] ** 2))

    reference_level = np.nanmedian(
        0.5 * (np.abs(uv_after[mask_after]) + np.abs(orange_after[mask_after]))
    )

    if np.isfinite(reference_level) and reference_level > 0:
        fractional_rms_before = rms_before / reference_level
        fractional_rms_after = rms_after / reference_level
    else:
        fractional_rms_before = np.nan
        fractional_rms_after = np.nan

    if np.isfinite(no_shift_fit["cost"]) and best_fit["cost"] > 0:
        improvement_factor = no_shift_fit["cost"] / best_fit["cost"]
    else:
        improvement_factor = np.nan

    shift_scan = np.linspace(shift_bounds[0], shift_bounds[1], scan_points)

    cost_scan = np.full(shift_scan.shape, np.nan, dtype=float)
    scale_scan = np.full(shift_scan.shape, np.nan, dtype=float)
    n_valid_scan = np.zeros(shift_scan.shape, dtype=int)

    for index, trial_shift in enumerate(shift_scan):
        trial = evaluate_delta(trial_shift)

        cost_scan[index] = trial["cost"]
        scale_scan[index] = trial["scale"]
        n_valid_scan[index] = np.sum(trial["valid"])

    return {
        "uv": uv_output,
        "orange": orange_output,
        "fit": {
            "delta_lambda": delta_lambda,
            "sky_orange_scale": orange_scale,
            "science_orange_scale": science_scale,
            "orange_scale": orange_scale,
            "orange_scale_no_shift": (orange_scale_no_shift),
            "cost": best_fit["cost"],
            "cost_no_shift": no_shift_fit["cost"],
            "improvement_factor": (improvement_factor),
            "rms_before": rms_before,
            "rms_after": rms_after,
            "fractional_rms_before": (fractional_rms_before),
            "fractional_rms_after": (fractional_rms_after),
            "n_fit": int(np.sum(mask_after)),
            "n_nominal": int(n_base),
            "fit_range": (overlap_min, overlap_max),
            "orange_fit_min": orange_fit_min,
            "repair_orange_response": repair_orange_response,
            "shift_bounds": shift_bounds,
            "throughput_floor": (throughput_floor),
            "success": bool(optimization.success),
            "message": optimization.message,
        },
        "diagnostic": {
            "wave": grid,
            "mask_before": mask_before,
            "mask_after": mask_after,
            "weight": weight,
            "uv_sky_before": uv_before,
            "orange_sky_before": orange_before,
            "uv_sky_after": uv_after,
            "orange_sky_after": orange_after,
            "residual_before": residual_before,
            "residual_after": residual_after,
            "uv_response": response_uv_nominal,
            "orange_response": (response_orange_nominal),
            "uv_throughput_proxy": (throughput_uv),
            "orange_throughput_proxy": (throughput_orange),
        },
        "scan": {
            "delta_lambda": shift_scan,
            "cost": cost_scan,
            "orange_scale": scale_scan,
            "sky_orange_scale": scale_scan,
            "n_valid": n_valid_scan,
        },
    }


def _plot_fit(result):
    """Plot the fit diagnostics and corrected science spectra."""
    fit = result["fit"]
    diag = result["diagnostic"]
    scan = result["scan"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    mask_before = diag["mask_before"]
    mask_after = diag["mask_after"]
    uv = result["uv"]
    orange = result["orange"]
    plot_min = fit["fit_range"][0] - 155.0
    plot_max = fit["fit_range"][1] + 155.0
    uv_mask = (uv["wave"] >= plot_min) & (uv["wave"] <= plot_max)
    orange_mask = (orange["wave"] >= plot_min) & (orange["wave"] <= plot_max)

    ax = axes[0, 0]
    ax.plot(
        diag["wave"][mask_before], diag["uv_sky_before"][mask_before], label="UV before"
    )
    ax.plot(
        diag["wave"][mask_before],
        diag["orange_sky_before"][mask_before],
        label="Orange before",
    )
    ax.plot(
        diag["wave"][mask_after],
        diag["uv_sky_after"][mask_after],
        "--",
        label="UV after",
    )
    ax.plot(
        diag["wave"][mask_after],
        diag["orange_sky_after"][mask_after],
        "--",
        label="Orange after",
    )
    ax.set_xlabel("Wavelength [Å]")
    ax.set_ylabel("Sky flux")
    ax.set_title("Sky overlap")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(
        uv["wave"][uv_mask],
        uv["sky_original"][uv_mask],
        alpha=0.5,
        label="UV sky before",
    )
    ax.plot(
        orange["wave"][orange_mask],
        orange["sky_original"][orange_mask],
        alpha=0.5,
        label="Orange sky before",
    )
    ax.plot(uv["wave"][uv_mask], uv["sky"][uv_mask], "--", label="UV sky after")
    ax.plot(
        orange["wave"][orange_mask],
        orange["sky"][orange_mask],
        "--",
        label="Orange sky after",
    )
    ax.set_xlabel("Wavelength [Å]")
    ax.set_ylabel("Sky flux")
    ax.set_title("Sky spectrum near dichroic")
    ax.legend()

    ax = axes[1, 0]
    good = np.isfinite(scan["cost"])
    ax.plot(scan["delta_lambda"][good], scan["cost"][good])
    ax.axvline(
        fit["delta_lambda"],
        ls="--",
        label=rf"Best $\Delta\lambda$ = {fit['delta_lambda']:+.3f} Å",
    )
    ax.axvline(0.0, ls=":", alpha=0.6)
    ax.set_xlabel(r"Dichroic shift $\Delta\lambda$ [Å]")
    ax.set_ylabel("Weighted mismatch")
    ax.set_title("Dichroic-shift objective")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(
        uv["wave"][uv_mask],
        uv["science_original"][uv_mask],
        alpha=0.5,
        label="UV science before",
    )
    ax.plot(
        orange["wave"][orange_mask],
        orange["science_original"][orange_mask],
        alpha=0.5,
        label="Orange science before",
    )
    ax.plot(uv["wave"][uv_mask], uv["science"][uv_mask], "--", label="UV science after")
    ax.plot(
        orange["wave"][orange_mask],
        orange["science"][orange_mask],
        "--",
        label="Orange science after",
    )
    ax.set_xlabel("Wavelength [Å]")
    ax.set_ylabel("Science flux")
    ax.set_title("Science spectrum near dichroic")
    ax.legend()

    title = (
        f"Dichroic shift = {fit['delta_lambda']:+.3f} Å   |   "
        f"Sky orange scale = {fit['sky_orange_scale']:.4f}   |   "
        f"Science orange scale = {fit['science_orange_scale']:.4f}   |   "
        f"Fractional RMS: {fit['fractional_rms_before']:.4f} → "
        f"{fit['fractional_rms_after']:.4f}"
    )
    fig.suptitle(title, fontsize=14)
    plt.savefig("dichroic_uv_orange_fit.png", dpi=150)


def _plot_response_diagnostic(result):
    """Plot nominal and shifted response functions for both channels."""
    fit = result["fit"]
    uv = result["uv"]
    orange = result["orange"]
    delta = fit["delta_lambda"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    plot_min = 4615.0
    plot_max = 4670.0

    def plot_channel_response(
        ax_response, ax_correction, channel, name, include_repaired=False
    ):
        """Plot response and correction diagnostics for one channel."""
        wave = channel["wave"]
        response_original = channel["response_original"]
        response_repaired = channel.get("response_repaired", response_original)
        response_shifted = channel["response_shifted"]

        response_correction = np.full_like(response_original, np.nan)
        np.divide(
            response_shifted,
            response_original,
            out=response_correction,
            where=response_original != 0,
        )

        sampled_wave = wave - delta
        extrapolated = (sampled_wave < np.nanmin(wave)) | (
            sampled_wave > np.nanmax(wave)
        )
        plot_mask = (
            (wave >= plot_min)
            & (wave <= plot_max)
            & np.isfinite(response_original)
            & np.isfinite(response_shifted)
        )

        reference_mask = (
            (wave >= 4640.0)
            & (wave <= 4645.0)
            & np.isfinite(response_original)
            & (response_original != 0)
        )
        if np.any(reference_mask):
            reference = np.nanmedian(response_original[reference_mask])
        else:
            reference = 1.0
        if not np.isfinite(reference) or reference == 0:
            reference = 1.0

        ax_response.plot(
            wave[plot_mask],
            response_original[plot_mask] / reference,
            label=f"{name} original",
        )
        if include_repaired:
            ax_response.plot(
                wave[plot_mask],
                response_repaired[plot_mask] / reference,
                "-.",
                label=f"{name} repaired",
            )
        ax_response.plot(
            wave[plot_mask],
            response_shifted[plot_mask] / reference,
            "--",
            label=f"{name} shifted",
        )

        extrapolated_plot = plot_mask & extrapolated
        if np.any(extrapolated_plot):
            ax_response.axvspan(
                np.nanmin(wave[extrapolated_plot]),
                np.nanmax(wave[extrapolated_plot]),
                alpha=0.15,
                label="Extrapolated",
            )

        correction_mask = plot_mask & np.isfinite(response_correction)
        ax_correction.plot(
            wave[correction_mask], response_correction[correction_mask], label=name
        )

        if np.any(extrapolated_plot):
            ax_correction.axvspan(
                np.nanmin(wave[extrapolated_plot]),
                np.nanmax(wave[extrapolated_plot]),
                alpha=0.15,
            )

    plot_channel_response(axes[0, 0], axes[1, 0], uv, "UV")
    axes[0, 0].set_title("UV response")
    axes[0, 0].set_xlabel("Wavelength [Å]")
    axes[0, 0].set_ylabel("Normalized response")
    axes[0, 0].legend()

    axes[1, 0].axhline(1.0, ls=":", alpha=0.6)
    axes[1, 0].set_title("UV response correction")
    axes[1, 0].set_xlabel("Wavelength [Å]")
    axes[1, 0].set_ylabel(r"$R(\lambda-\Delta\lambda) / R(\lambda)$")

    plot_channel_response(
        axes[0, 1], axes[1, 1], orange, "Orange", include_repaired=True
    )
    axes[0, 1].set_title("Orange response")
    axes[0, 1].set_xlabel("Wavelength [Å]")
    axes[0, 1].set_ylabel("Normalized response")
    axes[0, 1].legend()

    axes[1, 1].axhline(1.0, ls=":", alpha=0.6)
    axes[1, 1].set_title("Orange response correction")
    axes[1, 1].set_xlabel("Wavelength [Å]")
    axes[1, 1].set_ylabel(r"$R(\lambda-\Delta\lambda) / R(\lambda)$")

    fig.suptitle(
        f"Response-function effect for $\\Delta\\lambda={delta:+.3f}$ Å", fontsize=14
    )
    fig.savefig("dichroic_uv_orange_response.png", dpi=150)


def main():
    """Fit and plot the dichroic correction for two LRS2 spectra."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("uv_filename", help="UV spectrum filename")
    parser.add_argument("orange_filename", help="Orange spectrum filename")
    parser.add_argument(
        "--repair-orange-response",
        action="store_true",
        help="Repair the blue edge of the Orange response before fitting.",
    )
    args = parser.parse_args()

    result = fit_lrs2_dichroic(
        args.uv_filename,
        args.orange_filename,
        fit_range=(4635.0, 4645.0),
        # Orange data below 4638 Angstrom is excluded from the fit.
        orange_fit_min=4638.0,
        repair_orange_response=args.repair_orange_response,
        shift_bounds=(-10.0, 12.0),
        scale_bounds=(0.95, 1.05),
        throughput_floor=0.02,
    )
    _plot_fit(result)
    _plot_response_diagnostic(result)


if __name__ == "__main__":
    main()
