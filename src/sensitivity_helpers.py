import numpy as np
import pandas as pd

from model_config import alpha_rate, baseline_params, teredo_k
from decay_model import compute_total_decay_full, build_rate_scale


def scenario_total_and_rate(ds, material, use_teredo, steel_ref_ref, alpha=alpha_rate):
    return build_rate_scale(
        ds,
        material=material,
        use_teredo=use_teredo,
        steel_ref_ref=steel_ref_ref,
        alpha=alpha
    )


def summarise_ratio(
    test_rate,
    base_rate,
    parameter,
    value,
    scenario,
    metric="rate_scale"
):
    ratio = (test_rate / base_rate).values
    ratio = ratio[np.isfinite(ratio)]

    return {
        "parameter": parameter,
        "value": value,
        "scenario": scenario,
        "metric": metric,
        "median_ratio": np.nanmedian(ratio),
        "p05_ratio": np.nanpercentile(ratio, 5),
        "p25_ratio": np.nanpercentile(ratio, 25),
        "p75_ratio": np.nanpercentile(ratio, 75),
        "p95_ratio": np.nanpercentile(ratio, 95),
        "min_ratio": np.nanmin(ratio),
        "max_ratio": np.nanmax(ratio)
    }


def run_scenario(
    ABMsubset,
    steel_ref_ref,
    material,
    teredo="none",
    params=None
):
    if params is None:
        params = baseline_params.copy()

    k = teredo_k[teredo]
    use_teredo = (teredo != "none") and (material in ["wood", "mixed"])

    ds = compute_total_decay_full(
        ABMsubset.copy(),
        wreck_material=material,
        k_teredo=k,
        **params
    )

    CR_total, rate_scale = build_rate_scale(
        ds,
        material=material,
        use_teredo=use_teredo,
        steel_ref_ref=steel_ref_ref,
        alpha=alpha_rate
    )

    return ds, CR_total, rate_scale
