"""The DIII-D configuration file."""

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
            "range": [2, 7.5],
            "axis_name": r"""$q_{95}$""",
        },
        "ip": {
            "range": [-1.5e6, 1.5e6],
            "axis_name": "$I_p$",
        },
        "z_error": {
            "range": [-0.1, 0.1],
            "axis_name": r"""$z_{error}$""",
        },
        "n_over_ncrit": {
            "range": [-1, 2],
            "axis_name": r"""$n/n_{crit}$""",
        },
    },
}
