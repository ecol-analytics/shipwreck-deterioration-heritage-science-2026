alpha_rate = 0.7
support_age_max = 125.0
tail_lambda0 = 0.01
max_threshold_years = 300

baseline_params = dict(
    Ea=50000.0,
    U_ref=0.3,
    K_O2=30.0,
    S_mid_corr=12.0,
    S_steep_corr=6.0,
    burial_alpha=0.7,
    burial_strength=0.5
)

teredo_k = {
    "none": 0.0,
    "low": 0.02,
    "medium": 0.10,
    "high": 0.50
}

scenarios = [
    ("steel", "none"),
    ("wood", "none"),
    ("wood", "low"),
    ("wood", "medium"),
    ("wood", "high"),
    ("mixed", "none"),
    ("mixed", "low"),
    ("mixed", "medium"),
    ("mixed", "high"),
]

thresholds = ["2.0", "2.5", "3.0", "3.5", "3.9"]
