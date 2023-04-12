"""vis/disruptivity_vis.py contains all the plotting routines used to plot 1D and 2D disruptivity
maps."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from vis.plot_helpers import *

# TODO: Potentially fix the duplication of code between disruptivity1D and disruptivity2D when we introduce classes.
# TODO: When we get the big dictionaries and classes going, we should fix axis titles.

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

    assert "range" in entry, f"Entry {key} of entry_dict missing range field."
    ax.set_xlim(*entry["range"])

    # And the axis names
    assert (
        "axis_name" in entry
    ), f"Entry {key} of entry_dict missing axis_name field."
    ax.set_xlabel(entry["axis_name"])

    # Plot the bar plot
    width = (bins[0][1] - bins[0][0]) / 1.1
    ax.bar(
        bins[0][0:-1],
        disruptivity,
        color="firebrick",
        yerr=error,
        width=width,
        align="center",
        zorder=2,
        edgecolor=None,
        capsize=3,
    )

    ax.set_yscale("log")
    ax.grid(visible=None, which="both", axis="both")
    ax.set_ylabel("Disruptivity [1/s]")


# ----- 2D Disruptivity ----


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

    # Heatmap
    cax = ax.imshow(
        disruptivity.T,
        cmap="viridis",
        origin="lower",
        # interpolation="spline16",
        aspect="auto",
        extent=extent,
        norm=colors.LogNorm(),
    )

    # Axis Titles
    ax.set_xlabel(axis_name_list[0])
    ax.set_ylabel(axis_name_list[1])

    # Colorbar
    cbar = plt.colorbar(cax, label="Disruptivity ($s^{-1}$)")
    cbar.ax.tick_params(labelsize="large")
    cbar.set_label(label="Disruptivity ($s^{-1}$)", size="large")

    # TODO: Is there some way to visualize the errors in this type of plot?
