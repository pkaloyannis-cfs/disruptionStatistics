"""vis/disruptivity_vis.py contains all the plotting routines used to plot 1D and 2D disruptivity
maps."""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

figure_type = matplotlib.figure.Figure
subplot_type = matplotlib.axes._axes.Axes
figure_subplot_type = tuple[figure_type, subplot_type]

# TODO: Potentially fix the duplication of code between disruptivity1D and disruptivity2D when we introduce classes.
# TODO: When we get the big dictionaries and classes going, we should fix axis titles.

# ---- 1D Disruptivity ----


def plot_disruptivity1d(
    fpath: str, disruptivity: np.ndarray, error: np.ndarray, bins: np.ndarray
) -> figure_subplot_type:
    """Plots and saves a 1D disruptivity profile.

    Args:
        fpath (str): The filepath to save the plot to.
        disruptivity (np.ndarray): The disruptivity histogram of size (n_bins-1).
        error (np.ndarray): The error histogram of size (n_bins-1).
        bins (np.ndarray): The bin array of size (nbins)

    Returns:
        fig (plt.figure.Figure): Figure for plotting.
        ax (plt.axes._subplots.AxesSubplot): Axes for plotting.
    """

    # Prep the figure.
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Plot the figure.
    subplot_disruptivity1d(ax, disruptivity, error, bins)

    # Save and return the figure.
    fig.savefig(fpath, dpi=400, facecolor="w", bbox_inches="tight")
    return (fig, ax)


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


def plot_disruptivity2d(
    fpath: str,
    disruptivity: np.ndarray,
    error: np.ndarray,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
) -> figure_subplot_type:
    """Plots and saves a 2D disruptivity profile.

    Args:
        fpath (str): The filepath to save the plot to.
        disruptivity (np.ndarray): The disruptivity histogram of size (ny_bins-1, nx_bins-1).
        error (np.ndarray): The error histogram of size (ny_bins-1, nx_bins-1).
        x_bins (np.ndarray): The x bin array of size (nx_bins).
        y_bins (np.ndarray): The y bin array of size (ny_bins).

    Returns:
        fig (plt.figure.Figure): Figure for plotting.
        ax (plt.axes._subplots.AxesSubplot): Axes for plotting.
    """

    # Prep the figure.
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Plot the figure.
    subplot_disruptivity2d(ax, disruptivity, error, x_bins, y_bins)

    # Save and return the figure.
    fig.savefig(fpath, dpi=400, facecolor="w", bbox_inches="tight")
    return (fig, ax)


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
