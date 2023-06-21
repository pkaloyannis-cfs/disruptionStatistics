"""The DIII-D configuration file."""

# Structure.
CONFIG = {
    # Setup Information
    "name": "DIII-D",
    "filename_prefix": "d3d",
    "dataloc": "../data/d3d-db-220420.mat",
    "entry_dict": {
        "kappa": {
            "range": [1.0, 2.2],
            "axis_name": "$\kappa$",
        },
        "beta_n": {
            "range": [0, 4],
            "axis_name": r"""$\beta_n$""",
        },
        "li": {
            "range": [0.4, 3],
            "axis_name": "$l_i$",
        },
        "q95": {
            "range": [1, 10],
            "axis_name": r"""$q_{95}$""",
        },
        "ip": {
            "range": [-2e6, 2e6],
            "axis_name": "$I_p$",
        },
        "n_e": {"range": [0, 1.2e20], "axis_name": "$n_e$"},
        "p_rad": {"range": [0, 2e7], "axis_name": r"""$P_{Rad}$"""},
        "n1rms_normalized": {
            "range": [0, 2e-3],
            "axis_name": "Normalized 2/1\n Rotating Mode\n Amplitude",
        },
        "n_equal_1_normalized": {
            "range": [0, 1e-3],
            "axis_name": "Normalized 2/1\n Locked Mode\n Amplitude",
        },
        "radiated_fraction": {
            "range": [0, 3],
            "axis_name": "Radiated Fraction",
        },
    },
}
