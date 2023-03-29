"""vis/probability_vis.py hold code that is useful for generating disruption probability plots."""
# Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vis.plot_helpers import *

# ---- Disrupt Prob Over Time ----


def subplot_disrupt_rate_over_time(
    ax: subplot_type,
    dataframe: pd.core.frame.DataFrame,
    index_dict: dict,
    filter_size=300,
):
    """Subplot of the disruption rate over a TCV like shot number axis.

    This plot is inspired by figure of DeVries 2009.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        dataframe (pd.core.frame.DataFrame): The tokamak disruptivity dataframe.
        index_dict (dict): The index dictionary for each index type.
        filter_size (int, optional): Running average filter size. Defaults to 300.
    """

    # Label each shot number as disrupted or not.
    # 1. Find all the shots that have disrupted
    # 2. Find the total number of shots
    # 3. Mark disrupted shots with 1, non disrupted shots with 0 via intersection.
    flattop_disrupt_shot_number = dataframe.shot[
        index_dict["indices_disrupt_time_in_flattop"]
    ]
    shots = np.unique(dataframe.shot)
    shot_array_pre_avg = np.isin(shots, flattop_disrupt_shot_number).astype(
        np.int32
    )

    # Define the running average filter
    conv_filter = np.ones(filter_size) / filter_size

    # Define the TCV style shot number and axis limits.
    shot_num_TCV_style = np.arange(len(shots))
    xlim = (0, shot_num_TCV_style[-1])
    ylim = (0, 1)

    # Compute the running average
    shot_array = np.convolve(shot_array_pre_avg, conv_filter, mode="same")

    # Plot the running average and filter size
    ax.plot(shot_num_TCV_style, shot_array, color="k")
    ax.text(
        1000,
        0.865,
        f"Filter Size: {filter_size}",
        bbox=dict(
            boxstyle="round",
            ec=(1.0, 0.5, 0.5),
            fc=(1.0, 0.8, 0.8),
        ),
    )

    # Plot the filter warmup
    ax.fill_between(
        [xlim[0], filter_size],
        [ylim[-1], ylim[-1]],
        [ylim[0], ylim[0]],
        color="g",
        alpha=0.3,
    )
    ax.fill_between(
        [xlim[1] - filter_size, xlim[1]],
        [ylim[-1], ylim[-1]],
        [ylim[0], ylim[0]],
        color="g",
        alpha=0.3,
    )

    # Grid, xlims, ylims
    ax.grid()
    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    ax.set_ylabel("Average Disruption Rate")
    ax.set_xlabel("TCV Style Shot Number")


# ---- Intentional Disruptivity Rate ----


def subplot_intentional_rate_over_time(
    ax: subplot_type,
    dataframe: pd.core.frame.DataFrame,
    index_dict: dict,
):
    """AI is creating summary for subplot_intentional_rate_over_time.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        dataframe (pd.core.frame.DataFrame): The tokamak disruptivity dataframe.
        index_dict (dict): The index dictionary for each index type.
    """
