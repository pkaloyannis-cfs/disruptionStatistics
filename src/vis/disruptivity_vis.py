"""vis/disruptivity_vis.py contains all the plotting routines used to plot 1D and 2D disruptivity
maps."""

import numpy as np
import matplotlib.pyplot as plt
from vis.plot_helpers import *

# TODO: Potentially fix the duplication of code between disruptivity1D and disruptivity2D when we introduce classes.
# TODO: When we get the big dictionaries and classes going, we should fix axis titles.

# ---- 1D Disruptivity ----


def subplot_disruptivity1d(
    ax: subplot_type,
    disruptivity: np.ndarray,
    error: np.ndarray,
    bins: np.ndarray,
):
    """Creates a subplot of a 1D disruptivity map.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        disruptivity (np.ndarray): The disruptivity histogram of size (n_bins-1).
        error (np.ndarray): The error histogram of size (n_bins-1).
        bins (np.ndarray): The bin array of size (nbins).
    """

    # Asserts to ensure these are the correct bins and arrays
    assert len(disruptivity) + 1 == len(bins), (
        "Bin and Disruptivity array size mismatch. Expected bin size:"
        f" {len(disruptivity)+1}, Actual bin size: {len(bins)}"
    )
    assert len(disruptivity) == len(error), (
        "Error and Disruptivity array size mismatch. disruptivity size:"
        f" {len(disruptivity)}, error Size: {len(error)}"
    )

    # Plot the bar plot
    width = (bins[1] - bins[0]) / 1.1
    ax.bar(
        bins[0:-1],
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
    x_bins: np.ndarray,
    y_bins: np.ndarray,
):
    """Creates a subplot of a 2D disruptivity map.

    Args:
        ax (subplot_type): The matplotlib axes to plot on.
        disruptivity (np.ndarray): The disruptivity histogram of size (ny_bins-1, nx_bins-1).
        error (np.ndarray): The error histogram of size (ny_bins-1, nx_bins-1).
        x_bins (np.ndarray): The x bin array of size (nx_bins).
        y_bins (np.ndarray): The y bin array of size (ny_bins).
    """

    # Heatmap
    cax = ax.imshow(
        disruptivity.T,
        cmap="viridis",
        origin="lower",
        interpolation="spline16",
        aspect="auto",
        extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]],
    )

    # Colorbar
    cbar = plt.colorbar(cax, label="Disruptivity ($s^{-1}$)")
    cbar.ax.tick_params(labelsize="large")
    cbar.set_label(label="Disruptivity ($s^{-1}$)", size="large")

    # TODO: Is there some way to visualize the errors in this type of plot?
