"""Shared linear-wave calculations used by Navigator and the MSM pipeline.

This module contains physical calculations only. It deliberately does not load
files, interpolate grids, calculate annual descriptors, or apply decay-model
response functions.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

GRAVITY_M_S2 = 9.81


def solve_wavenumber(
    wave_period_s: np.ndarray,
    water_depth_m: np.ndarray,
    *,
    gravity_m_s2: float = GRAVITY_M_S2,
    iterations: int = 20,
) -> np.ndarray:
    """Solve the linear-wave dispersion relation by Newton iteration.

    Solves ``omega**2 = g * k * tanh(k * h)`` for wavenumber ``k``.
    Invalid, non-positive period or depth values return NaN.
    """
    period, depth = np.broadcast_arrays(
        np.asarray(wave_period_s, dtype=float),
        np.asarray(water_depth_m, dtype=float),
    )
    valid = (
        np.isfinite(period)
        & np.isfinite(depth)
        & (period > 0.0)
        & (depth > 0.0)
    )

    angular_frequency = np.full_like(period, np.nan, dtype=float)
    angular_frequency[valid] = 2.0 * np.pi / period[valid]

    wavenumber = np.full_like(period, np.nan, dtype=float)
    wavenumber[valid] = angular_frequency[valid] ** 2 / gravity_m_s2

    for _ in range(iterations):
        kh = wavenumber[valid] * depth[valid]
        tanh_kh = np.tanh(kh)
        # cosh overflows harmlessly for very large kh, but clipping avoids the
        # warning and the derivative contribution is already effectively zero.
        kh_for_cosh = np.clip(kh, 0.0, 350.0)
        sech2_kh = 1.0 / np.cosh(kh_for_cosh) ** 2

        residual = (
            gravity_m_s2 * wavenumber[valid] * tanh_kh
            - angular_frequency[valid] ** 2
        )
        derivative = (
            gravity_m_s2 * tanh_kh
            + gravity_m_s2
            * wavenumber[valid]
            * depth[valid]
            * sech2_kh
        )
        wavenumber[valid] -= residual / derivative

    return wavenumber


def seabed_orbital_velocity_numpy(
    significant_wave_height_m: np.ndarray,
    mean_wave_period_s: np.ndarray,
    water_depth_m: np.ndarray,
    *,
    gravity_m_s2: float = GRAVITY_M_S2,
    iterations: int = 20,
) -> np.ndarray:
    """Return horizontal wave orbital-velocity amplitude at the seabed.

    Uses linear wave theory and the full finite-depth dispersion relation.
    Output units are m s-1. Invalid inputs remain NaN.
    """
    height, period, depth = np.broadcast_arrays(
        np.asarray(significant_wave_height_m, dtype=float),
        np.asarray(mean_wave_period_s, dtype=float),
        np.asarray(water_depth_m, dtype=float),
    )
    wavenumber = solve_wavenumber(
        period,
        depth,
        gravity_m_s2=gravity_m_s2,
        iterations=iterations,
    )

    valid = (
        np.isfinite(height)
        & np.isfinite(period)
        & np.isfinite(depth)
        & np.isfinite(wavenumber)
        & (height >= 0.0)
        & (period > 0.0)
        & (depth > 0.0)
    )

    velocity = np.full_like(height, np.nan, dtype=float)
    kh = np.clip(wavenumber[valid] * depth[valid], 0.0, 50.0)
    denominator = period[valid] * np.sinh(kh)
    velocity[valid] = np.pi * height[valid] / denominator
    return velocity


def seabed_orbital_velocity_xarray(
    significant_wave_height: xr.DataArray,
    mean_wave_period: xr.DataArray,
    water_depth: xr.DataArray,
    *,
    gravity_m_s2: float = GRAVITY_M_S2,
    iterations: int = 20,
) -> xr.DataArray:
    """Xarray/Dask wrapper for :func:`seabed_orbital_velocity_numpy`."""
    result = xr.apply_ufunc(
        seabed_orbital_velocity_numpy,
        significant_wave_height,
        mean_wave_period,
        water_depth,
        kwargs={
            "gravity_m_s2": gravity_m_s2,
            "iterations": iterations,
        },
        vectorize=False,
        dask="allowed",
        output_dtypes=[float],
    )
    result.name = "seabed_orbital_velocity"
    result.attrs.update(
        {
            "long_name": "wave-induced horizontal orbital velocity at the seabed",
            "units": "m s-1",
            "method": "linear wave theory with finite-depth dispersion relation",
        }
    )
    return result


# Temporary compatibility aliases for existing notebooks/scripts.
bottom_orbital_velocity_numpy = seabed_orbital_velocity_numpy
bottom_orbital_velocity_xr = seabed_orbital_velocity_xarray
