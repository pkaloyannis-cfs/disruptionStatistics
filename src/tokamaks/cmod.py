"""The DIII-D configuration file."""

# Structure.
CONFIG = {
    # Setup Information
    "name": "C-Mod",
    "filename_prefix": "cmod",
    "disrupt_warning_time_ms": 50,
    "disrupt_warning_window_ms": 10,
    "dataloc": "../data/CMod_disruption_warning_db.mat",
    "entry_dict": {
        "kappa": {
            "range": [0.8, 2.0],
            "axis_name": "$\kappa$",
        },
        "beta_n": {
            "range": [0, 1.5],
            "axis_name": "$\beta_n$",
        },
        "li": {
            "range": [0.75, 2.25],
            "axis_name": "$l_i$",
        },
        "q95": {
            "range": [2, 7.5],
            "axis_name": "$\beta_n$",
        },
        "ip": {
            "range": [-1.5e6, 1.5e6],
            "axis_name": "$I_p$",
        },
    },
}
