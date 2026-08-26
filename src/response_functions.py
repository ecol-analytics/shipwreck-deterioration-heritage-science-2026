import numpy as np


def salinity_corrosion_modifier(
    S,
    S_mid=12.0,
    S_steep=6.0,
    f_min=0.05,
    f_max=1.0
):
    raw = 1.0 / (1.0 + np.exp(-(S - S_mid) / S_steep))
    fS = f_min + (f_max - f_min) * raw
    return fS.clip(min=f_min, max=f_max)


def teredo_salinity_modifier(
    S,
    S_mid=9.0,
    S_steep=2.0,
    f_min=0.0,
    f_max=1.0
):
    raw = 1.0 / (1.0 + np.exp(-(S - S_mid) / S_steep))
    fS = f_min + (f_max - f_min) * raw
    return fS.clip(min=f_min, max=f_max)
