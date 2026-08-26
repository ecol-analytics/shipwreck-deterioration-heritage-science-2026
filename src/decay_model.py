import numpy as np
import xarray as xr

from wave_physics import seabed_orbital_velocity_xarray

from response_functions import (
    salinity_corrosion_modifier,
    teredo_salinity_modifier
)

from model_config import alpha_rate, support_age_max, tail_lambda0


def compute_total_decay_full(
    ds,
    wreck_material="steel",
    k_corr=0.020,
    k_erosion=0.020,
    k_bio=0.050,
    k_bio_steel=0.005,
    k_bio_wood=0.020,
    k_teredo=0.5,
    T_ref_C=20.0,
    Ea=50000.0,
    U_ref=0.3,
    U_current_ref=None,
    U_wave_ref=None,
    hydro_interaction_strength=0.0,
    K_O2=30.0,
    S_mid_corr=12.0,
    S_steep_corr=6.0,
    burial_alpha=0.7,
    burial_strength=0.5,
    steel_fraction_mixed=0.6,
    wood_fraction_mixed=0.4
):
    out = ds.copy()

    T_C = out["bottomT_mean"] if "bottomT_mean" in out else out["thetao_mean"]
    T_K = T_C + 273.15
    T_ref_K = T_ref_C + 273.15

    O2 = out["o2_mean"]
    PO4 = out["po4_mean"]
    NO3 = out["no3_mean"]
    NPPV = out["nppv_mean"]
    S = out["so_mean"]
    pH = out["ph_mean"]

    Fe = out["fe_mean"] if "fe_mean" in out else xr.full_like(T_C, 0.05)
    Chl = out["chl_mean"] if "chl_mean" in out else xr.full_like(T_C, 0.5)
    substrate = out["substrate_soft_hard"]

    # ----------------------------------------------------------
    # Current descriptors
    # ----------------------------------------------------------

    U_mean = out["uo_mean"]
    V_mean = out["vo_mean"]

    # Magnitude of the long-term mean vector:
    # useful for prevailing direction and directional persistence,
    # but not as the scalar current forcing.
    current_vector_speed = np.hypot(U_mean, V_mean)

    # Mean of monthly current-speed magnitudes:
    # appropriate scalar measure of hydrodynamic forcing.
    if "current_speed_mean" in out:
        current_speed = out["current_speed_mean"]
    else:
        # Backward-compatible fallback for older datasets
        current_speed = current_vector_speed

    # Mathematical direction of the mean vector:
    # degrees anticlockwise from east, in the range 0–360.
    current_direction_math = (
        np.degrees(np.arctan2(V_mean, U_mean)) + 360.0
    ) % 360.0

    # Compass bearing towards which the current flows:
    # 0° = north, 90° = east.
    current_direction_to = (
        90.0 - current_direction_math
    ) % 360.0

    # Directional persistence:
    # 1 = consistently aligned flow;
    # 0 = energetic but cancelling/reversing flow.
    eps = 1e-12
    directional_persistence = xr.where(
        current_speed > eps,
        current_vector_speed / current_speed,
        np.nan
    ).clip(0.0, 1.0)

    # ----------------------------------------------------------
    # Wave descriptors
    # ----------------------------------------------------------

    if "seabed_orbital_velocity_mean" in out:
        # Preferred production representation: seabed orbital velocity
        # calculated from the original 3-hourly wave fields before temporal
        # averaging. This avoids the nonlinear aggregation bias introduced by
        # calculating orbital velocity from long-term mean wave height and
        # period.
        orbital_velocity = out["seabed_orbital_velocity_mean"]
    else:
        # Backward-compatible fallback for older regional datasets that do not
        # yet contain the precomputed 3-hourly-derived orbital-velocity field.
        depth = out["deptho_mean"]
        wave_period = out["VTM02_mean_mean"].where(
            out["VTM02_mean_mean"] > 1e-6
        )
        wave_height = out["VHM0_mean_mean"].where(
            out["VHM0_mean_mean"] >= 0.0
        )

        orbital_velocity = seabed_orbital_velocity_xarray(
            significant_wave_height=wave_height,
            mean_wave_period=wave_period,
            water_depth=depth,
        )

    # Retained as a diagnostic for comparison with the previous formulation.
    # It no longer drives the hydrodynamic modifier directly.
    total_velocity = np.sqrt(
        current_speed**2 + orbital_velocity**2
    )

    # ----------------------------------------------------------
    # Regional hydrodynamic exposure
    # ----------------------------------------------------------
    #
    # U_ref is retained for backward compatibility with existing parameter
    # dictionaries and sensitivity scripts. Unless separate reference scales
    # are supplied, both current and wave terms use U_ref.

    if U_current_ref is None:
        U_current_ref = U_ref

    if U_wave_ref is None:
        U_wave_ref = U_ref

    if U_current_ref <= 0:
        raise ValueError("U_current_ref must be greater than zero.")

    if U_wave_ref <= 0:
        raise ValueError("U_wave_ref must be greater than zero.")

    if hydro_interaction_strength < 0:
        raise ValueError(
            "hydro_interaction_strength must be non-negative."
        )

    fCurrent = (
        current_speed
        / (current_speed + U_current_ref)
    ).clip(0.0, 1.0)

    fWave = (
        orbital_velocity
        / (orbital_velocity + U_wave_ref)
    ).clip(0.0, 1.0)

    # Independent contribution from either current or wave action.
    #
    # This bounded union-style formulation retains the full effect of one
    # component when the other is negligible, while avoiding simple
    # double-counting where both are elevated:
    #
    #     fHydro_base = fCurrent + fWave - fCurrent * fWave
    #
    # Strong current with negligible waves therefore remains strongly exposed,
    # which is important in deep current-dominated regions such as the
    # Norwegian Trench.

    fHydro_base = (
        fCurrent
        + fWave
        - fCurrent * fWave
    ).clip(0.0, 1.0)

    # Optional non-additive enhancement where both current and wave exposure
    # are elevated. A value of zero recovers the independent union model.
    hydro_compound_overlap = fCurrent * fWave

    fHydro = (
        fHydro_base
        * (
            1.0
            + hydro_interaction_strength
            * hydro_compound_overlap
        )
    ).clip(0.05, 1.0)

    R = 8.314
    fT = np.exp((-Ea / R) * ((1.0 / T_K) - (1.0 / T_ref_K))).clip(0.05, 5.0)
    fO2 = (O2 / (O2 + K_O2)).clip(0.05, 1.0)

    K_PO4, K_NO3, K_NPPV = 1.0, 5.0, 0.1
    fNUT = (
        (PO4 / (PO4 + K_PO4)) *
        (NO3 / (NO3 + K_NO3)) *
        (NPPV / (NPPV + K_NPPV))
    ).clip(0.05, 1.0)

    sigma_pH = 0.5
    fPH = np.exp(-0.5 * ((pH - 8.0) / sigma_pH) ** 2).clip(0.05, 1.0)

    fS_corr = salinity_corrosion_modifier(
        S,
        S_mid=S_mid_corr,
        S_steep=S_steep_corr,
        f_min=0.05,
        f_max=1.0
    )

    fFe = (Fe / (Fe + 0.1)).clip(0.05, 1.0)
    fChl = (Chl / (Chl + 1.0)).clip(0.05, 1.0)

    hardness = ((substrate - 1.0) / 9.0).clip(0.0, 1.0)
    fBurial = ((1.0 - hardness) ** burial_alpha).clip(0.05, 0.95)
    burial_modifier = (1.0 - burial_strength * fBurial).clip(0.05, 1.0)
    exposure = (fHydro * burial_modifier).clip(0.05, 1.0)

    out["fCurrent"] = fCurrent
    out["fWave"] = fWave
    out["fHydro_base"] = fHydro_base
    out["hydro_compound_overlap"] = hydro_compound_overlap
    out["fHydro"] = fHydro
    out["hydrodynamic_exposure"] = fHydro
    out["fBurial"] = fBurial
    out["burial_modifier"] = burial_modifier
    out["exposure"] = exposure
    out["fT"] = fT
    out["fO2"] = fO2
    out["fS_corr"] = fS_corr
    out["orbital_velocity"] = orbital_velocity
    out["total_velocity"] = total_velocity
    out["fNUT"] = fNUT
    out["fFe"] = fFe
    out["fChl"] = fChl
    out["fPH"] = fPH
    out["ph_mean"] = pH
    out["current_speed"] = current_speed
    out["current_vector_speed"] = current_vector_speed
    out["current_direction_to"] = current_direction_to
    out["directional_persistence"] = directional_persistence

    # Canonical descriptive names. Legacy names remain during migration so the
    # dashboard and sensitivity notebooks continue to run unchanged.
    out["mean_current_speed"] = current_speed
    out["mean_vector_current_speed"] = current_vector_speed
    out["current_directional_persistence"] = directional_persistence
    out["seabed_orbital_velocity"] = orbital_velocity
    out["current_forcing_modifier"] = fCurrent
    out["wave_forcing_modifier"] = fWave
    out["hydrodynamic_forcing_base"] = fHydro_base
    out["current_wave_overlap"] = hydro_compound_overlap
    out["hydrodynamic_forcing_modifier"] = fHydro
    out["seabed_exposure_modifier"] = exposure
    out["temperature_modifier"] = fT
    out["oxygen_modifier"] = fO2
    out["salinity_corrosion_modifier"] = fS_corr
    out["nutrient_modifier"] = fNUT
    out["iron_modifier"] = fFe
    out["chlorophyll_modifier"] = fChl
    out["ph_modifier"] = fPH

    CR_steel = k_corr * fT * fO2 * fS_corr * fHydro
    CR_wood_phys = k_erosion * exposure

    out["CR_steel_component"] = CR_steel.clip(min=1e-12)
    out["CR_wood_phys_component"] = CR_wood_phys.clip(min=1e-12)

    if wreck_material == "steel":
        CR_physics = CR_steel
    elif wreck_material == "wood":
        CR_physics = CR_wood_phys
    elif wreck_material == "mixed":
        CR_physics = steel_fraction_mixed * CR_steel + wood_fraction_mixed * CR_wood_phys
    else:
        CR_physics = CR_steel

    out["CR_physics"] = CR_physics.clip(min=1e-12)

    CR_bio_steel = (
        k_bio_steel * fT * fO2 * fNUT * fFe * fChl * exposure
    )

    CR_bio_wood = (
        k_bio_wood * fT * fO2 * exposure
    )

    CR_bio_steel = CR_bio_steel * (1.0 - 0.25 * fBurial)
    CR_bio_wood = CR_bio_wood * (1.0 - 0.50 * fBurial)

    if wreck_material == "steel":
        CR_bio = CR_bio_steel
    elif wreck_material == "wood":
        CR_bio = CR_bio_wood
    elif wreck_material == "mixed":
        CR_bio = steel_fraction_mixed * CR_bio_steel + wood_fraction_mixed * CR_bio_wood
    else:
        CR_bio = CR_bio_steel

    out["CR_bio"] = CR_bio.clip(min=1e-12)
    out["CR_bio_steel_component"] = CR_bio_steel.clip(min=1e-12)
    out["CR_bio_wood_component"] = CR_bio_wood.clip(min=1e-12)

    fS_T = teredo_salinity_modifier(
        S,
        S_mid=9.0,
        S_steep=2.0,
        f_min=0.0,
        f_max=1.0
    )

    T_opt, sigma_T = 20.0, 6.0
    fT_T = np.exp(-0.5 * ((T_C - T_opt) / sigma_T) ** 2).clip(0.0, 1.0)

    CR_teredo = k_teredo * fS_T * fT_T * fO2 * exposure

    if wreck_material == "steel":
        CR_teredo = xr.zeros_like(T_C)

    out["CR_teredo"] = CR_teredo.clip(min=0.0)

    out["CR_total_no_teredo"] = out["CR_physics"] + out["CR_bio"]
    out["CR_total_with_teredo"] = out["CR_total_no_teredo"] + out["CR_teredo"]
    out["CR_teredo_effect"] = out["CR_teredo"]

    out["steel_corrosion_rate"] = out["CR_steel_component"]
    out["wood_physical_decay_rate"] = out["CR_wood_phys_component"]
    out["physical_decay_rate"] = out["CR_physics"]
    out["biological_decay_rate"] = out["CR_bio"]
    out["teredo_decay_rate"] = out["CR_teredo"]
    out["total_decay_rate_without_teredo"] = out["CR_total_no_teredo"]
    out["total_decay_rate_with_teredo"] = out["CR_total_with_teredo"]

    out["log_CR_physics"] = np.log10(out["CR_physics"] + 1e-12)
    out["log_CR_bio"] = np.log10(out["CR_bio"] + 1e-12)

    out["fS_T"] = fS_T
    out["fT_T"] = fT_T

    return out


def teredo_to_k(level: str):
    from model_config import teredo_k

    if level not in teredo_k:
        return False, 0.0

    return level != "none", teredo_k[level]


def build_rate_scale(
    ds,
    material: str,
    use_teredo: bool,
    steel_ref_ref: float,
    alpha: float = alpha_rate
):
    CR_total = ds["CR_physics"] + ds["CR_bio"]

    if use_teredo and material in ["wood", "mixed"]:
        CR_total = CR_total + ds["CR_teredo"]

    rate_ratio = CR_total / (steel_ref_ref + 1e-12)
    rate_scale = (alpha * np.log1p(rate_ratio)).clip(min=0.05)

    return CR_total, rate_scale


def compute_steel_reference(ds, baseline_params):
    steel_ref = compute_total_decay_full(
        ds.copy(),
        wreck_material="steel",
        k_teredo=0.0,
        **baseline_params
    )

    steel_ref_rate = steel_ref["CR_physics"] + steel_ref["CR_bio"]
    steel_ref_ref = float(np.nanmedian(steel_ref_rate.values))

    return steel_ref, steel_ref_ref


def time_to_threshold_map(
    ds,
    material,
    use_teredo,
    steel_ref_ref,
    baseline_class,
    threshold=3.5,
    max_years=300,
    alpha: float = alpha_rate
):
    _, rate_scale = build_rate_scale(
        ds,
        material,
        use_teredo,
        steel_ref_ref,
        alpha=alpha
    )

    eff_time = np.zeros_like(rate_scale.values, dtype=float)
    class_field = np.full_like(rate_scale.values, np.nan, dtype=float)
    years_to_threshold = np.full_like(rate_scale.values, np.nan, dtype=float)

    for t in range(max_years + 1):
        eff_time += rate_scale.values

        in_empirical = eff_time <= support_age_max
        if np.any(in_empirical):
            class_field[in_empirical] = baseline_class(eff_time[in_empirical])

        in_tail = eff_time > support_age_max
        if np.any(in_tail):
            newly_tail = in_tail & np.isnan(class_field)
            if np.any(newly_tail):
                class_field[newly_tail] = baseline_class(support_age_max)

            k_eff = tail_lambda0 * rate_scale.values[in_tail]
            class_field[in_tail] = 4.0 - (4.0 - class_field[in_tail]) * np.exp(-k_eff)

        class_field = np.clip(class_field, 1.0, 4.0)

        newly_reached = (class_field >= threshold) & np.isnan(years_to_threshold)
        years_to_threshold[newly_reached] = t

    return xr.DataArray(
        years_to_threshold,
        coords={"latitude": ds["latitude"], "longitude": ds["longitude"]},
        dims=("latitude", "longitude"),
        name=f"time_to_{threshold}"
    )
