"""The C-MOD configuration file."""

# Structure.
CONFIG = {
    # Setup Information
    "name": "C-Mod",
    "filename_prefix": "cmod",
    "dataloc": "../data/CMod_disruption_warning_db.mat",
    "entry_dict": {
        "kappa": {
            "range": [0.8, 2.0],
            "axis_name": "$\kappa$",
        },
        "beta_n": {
            "range": [0, 1.5],
            "axis_name": r"""$\beta_n$""",
        },
        "li": {
            "range": [0.75, 2.25],
            "axis_name": "$l_i$",
        },
        "q95": {
            "range": [1, 7.5],
            "axis_name": r"""$q_{95}$""",
        },
        "ip": {
            "range": [-1.5e6, 1.5e6],
            "axis_name": "$I_p$ (MA)",
        },
        "z_error": {
            "range": [-0.1, 0.1],
            "axis_name": r"""$z_{error}$""",
        },
        "n_over_ncrit": {
            "range": [0, 2],
            "axis_name": r"""$n/n_{crit}$""",
        },
        "n_e": {"range": [0, 5e20], "axis_name": "$n_e$"},
        "p_rad": {"range": [0, 1e7], "axis_name": r"""$P_{Rad}$"""},
        "Mirnov_norm_btor": {
            "range": [0, 8],
            "axis_name": "Normalized Rotating $n=1$ Amplitude",
        },
        "n_equal_1_normalized": {
            "range": [0, 0.005],
            "axis_name": "Normalized Locked $n=1$ Amplitude",
        },
        "radiated_fraction": {
            "range": [0, 3],
            "axis_name": "Radiated Fraction",
        },
    },
}
