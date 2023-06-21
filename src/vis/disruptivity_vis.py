"""# Description

This file contains plotting routines related to disruptivity maps.

Note: No function in here is unit tested.

# Functions
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local
from vis.plot_helpers import *

# TODO: Potentially fix the duplication of code between disruptivity1D and disruptivity2D when we introduce classes.
# TODO: Make size asserts more informative for the plotters and entry dict lengths.

# ---- 1D Disruptivity ----


def subplot_disruptivity1d(
    ax: subplot_type,
    disruptivity: np.ndarray,
    error: np.ndarray,
    bins: np.ndarray,
    entry_dict: dict,
):
    """Creates a subplot of a 1D disruptivity map.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        disruptivity (np.ndarray): The disruptivity histogram of size (n_bins-1).
        error (np.ndarray): The error histogram of size (n_bins-1).
        bins (list): The list of bin arrays of size (nbins).
        entry_dict (list): List of the entries being plotted.
    """

    # Asserts to ensure these are the correct bins and arrays
    assert len(disruptivity) + 1 == len(bins[0]), (
        "Bin and Disruptivity array size mismatch. Expected bin size:"
        f" {len(disruptivity)+1}, Actual bin size: {len(bins[0])}"
    )
    assert len(disruptivity) == len(error), (
        "Error and Disruptivity array size mismatch. disruptivity size:"
        f" {len(disruptivity)}, error Size: {len(error)}"
    )
    assert len(entry_dict) == 1, "Too many entry_dict dimensions."

    # Get the limits and the axis name
    entry = entry_dict[list(entry_dict)[0]]

    # And the axis names
    assert (
        "axis_name" in entry
    ), f"Entry {key} of entry_dict missing axis_name field."
    ax.set_xlabel(entry["axis_name"])

    # Plot the bar plot
    width = bins[0][1] - bins[0][0]
    ax.bar(
        bins[0][0:-1] + width / 2,
        disruptivity,
        color="firebrick",
        yerr=error,
        width=width / 1.1,
        align="center",
        zorder=2,
        edgecolor=None,
        capsize=3,
    )

    ax.set_yscale("log")
    ax.grid(visible=None, which="both", axis="both")
    ax.set_ylabel("Disruptivity [1/s]")


# ----- 2D Disruptivity ----


# TODO: DOCSTRING HERE
def subplot_draw_trajectory(ax, dataframe, entry_dict, indices, shot):
    """AI is creating summary for subplot_draw_trajectory.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        dataframe (pd.DataFrame): The tokamak data dataframe.
        entry_dict (dict): The entry dictionary. Must have 2 keys.
        indices (np.ndarray): The list of indices of note.
        shot (int): The shot number to plot/

    Returns:
        (float): The start time of the flattop.
    """

    # Filter the datafram data to get the pulse of interest in the flattop
    shotlist_bool = np.isin(dataframe.shot, [shot])
    shot_indices = dataframe[shotlist_bool].index
    overlap = np.array(np.intersect1d(indices, shot_indices))
    print(f"Flat Top Starts: {dataframe.time[overlap[0]]} s")

    # Get the entry information
    entries = list(entry_dict.keys())
    assert len(entries) == 2, (
        "Incorrect number of entries for plotting trajectory over"
        + "2D disruptivity plot. Requires exactly 2 entries."
    )

    # Create line segments to color individually
    # Based on this example: https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
    time = dataframe["time"][overlap]
    x = dataframe[entries[0]][overlap]
    y = dataframe[entries[1]][overlap]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(time.min(), time.max())
    lc = LineCollection(
        segments, cmap="cool", norm=norm, capstyle="round"
    )  # solid_capstyle='round'
    # Set the values used for colormapping
    lc.set_array(time)
    lc.set_linewidth(2)

    # Add the collection
    line = ax.add_collection(lc)

    # RHS colorbar for time
    cbar = plt.colorbar(line, ax=ax, location="left", pad=0.13)
    cbar.ax.tick_params()
    cbar.set_label(label="Pulse Time (s)")

    # Re-apply axis labels and limits if needed.
    ax.set_xlim(entry_dict[entries[0]]["range"])
    ax.set_ylim(entry_dict[entries[1]]["range"])
    ax.set_xlabel(entry_dict[entries[0]]["axis_name"])
    ax.set_ylabel(entry_dict[entries[1]]["axis_name"])

    return dataframe.time[overlap[0]]


def subplot_grad2d(
    ax,
    disruptivity: np.ndarray,
    error: np.ndarray,
    bins: np.ndarray,
    entry_dict: dict,
):
    """Creates a subplot of a 2D disruptivity map.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        disruptivity (np.ndarray): The disruptivity histogram of size (ny_bins-1, nx_bins-1).
        error (np.ndarray): The error histogram of size (ny_bins-1, nx_bins-1).
        bins (list): List of the bin edges from the histogram.
        entry_dict (list): List of the entries being plotted.
    """

    # Asserts
    assert len(entry_dict) == 2, "Too many entry_dict dimensions."

    # Parse the dict
    extent = []
    axis_name_list = []
    for key in entry_dict:
        entry = entry_dict[key]

        # Follow the order of the dictionary to find x and y
        # Make sure the range is there
        assert (
            "range" in entry
        ), f"Entry {key} of entry_dict missing range field."
        extent.extend(entry["range"])

        # And the axis names
        assert (
            "axis_name" in entry
        ), f"Entry {key} of entry_dict missing axis_name field."
        axis_name_list.append(entry["axis_name"])

    # The Normal Heatmap
    cax = ax.imshow(
        disruptivity.T,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=extent,
    )

    # Axis Titles
    ax.set_xlabel(axis_name_list[0])
    ax.set_ylabel(axis_name_list[1])

    # Colorbar
    cbar = plt.colorbar(cax)
    cbar.ax.tick_params(labelsize="large")
    cbar.set_label(label="Disruptivity Derivative", size="large")


def subplot_disruptivity2d(
    ax: subplot_type,
    disruptivity: np.ndarray,
    error: np.ndarray,
    bins: np.ndarray,
    entry_dict: dict,
):
    """Creates a subplot of a 2D disruptivity map.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        disruptivity (np.ndarray): The disruptivity histogram of size (ny_bins-1, nx_bins-1).
        error (np.ndarray): The error histogram of size (ny_bins-1, nx_bins-1).
        bins (list): List of the bin edges from the histogram.
        entry_dict (list): List of the entries being plotted.
    """

    # Asserts
    assert len(entry_dict) == 2, "Too many entry_dict dimensions."

    # Parse the dict
    extent = []
    axis_name_list = []
    for key in entry_dict:
        entry = entry_dict[key]

        # Follow the order of the dictionary to find x and y
        # Make sure the range is there
        assert (
            "range" in entry
        ), f"Entry {key} of entry_dict missing range field."
        extent.extend(entry["range"])

        # And the axis names
        assert (
            "axis_name" in entry
        ), f"Entry {key} of entry_dict missing axis_name field."
        axis_name_list.append(entry["axis_name"])

    # The Normal Heatmap
    cax = ax.imshow(
        disruptivity.T,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=extent,
        norm=colors.LogNorm(),
    )

    # Now the masked values
    no_disruptions = np.ma.masked_where(
        disruptivity.T != -2, np.ones(disruptivity.T.shape)
    )
    all_disruptions = np.ma.masked_where(
        disruptivity.T != -3, np.ones(disruptivity.T.shape)
    )

    # The masked values
    # When no disruptions, draw a black box
    ax.imshow(
        no_disruptions,
        cmap="Accent",
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=0,
        vmax=1,
    )

    # When all disruptions, draw a yellow box
    ax.imshow(
        all_disruptions,
        cmap="Set1",
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=0,
        vmax=1,
    )

    # Axis Titles
    ax.set_xlabel(axis_name_list[0])
    ax.set_ylabel(axis_name_list[1])

    # Colorbar
    cbar = plt.colorbar(cax)
    cbar.ax.tick_params(labelsize="large")
    cbar.set_label(label="Disruptivity ($s^{-1}$)")

    # TODO: Is there some way to visualize the errors in this type of plot?


# Data Processing Visualizations
def plot_data_selection(
    disruptivity: np.ndarray,
    error: np.ndarray,
    bins: np.ndarray,
    entry_dict: dict,
):
    """Plots the results of the data selection routines.

    Args:
        disruptivity (np.ndarray): The disruptivity histogram of size (ny_bins-1, nx_bins-1).
        error (np.ndarray): The error histogram of size (ny_bins-1, nx_bins-1).
        bins (list): List of the bin edges from the histogram.
        entry_dict (list): List of the entries being plotted.
    """

    # Create the plot
    fig, ax = plt.subplots(1, 4, figsize=(4 * 4, 3))

    # Parse the dict
    extent = []
    axis_name_list = []
    for key in entry_dict:
        entry = entry_dict[key]

        # Follow the order of the dictionary to find x and y
        # Make sure the range is there
        assert (
            "range" in entry
        ), f"Entry {key} of entry_dict missing range field."
        extent.extend(entry["range"])

        # And the axis names
        assert (
            "axis_name" in entry
        ), f"Entry {key} of entry_dict missing axis_name field."
        axis_name_list.append(entry["axis_name"])

    title_list = ["Full Data", "No Data", "None Disrupted", "All Disrupted"]
    ans_list = [
        disruptivity >= 0,
        disruptivity == -1,
        disruptivity == -2,
        disruptivity == -3,
    ]
    for i, axis in enumerate(ax):
        axis.imshow(
            ans_list[i].T,
            cmap="viridis",
            origin="lower",
            aspect="auto",
            extent=extent,
        )
        axis.set_title(title_list[i])

        # Axis Titles
        axis.set_xlabel(axis_name_list[0])
        axis.set_ylabel(axis_name_list[1])

    return fig, ax
