"""# Description

This file contains general plotting helpers.

# Functions
"""

import matplotlib
import matplotlib.pyplot as plt
from collections.abc import Callable

# Type Definitions
figure_type = matplotlib.figure.Figure
subplot_type = matplotlib.axes._axes.Axes
figure_subplot_type = tuple[figure_type, subplot_type]


def plot_subplot(
    fpath: str, subplot_func: Callable, args: list, figsize=(4, 3)
) -> figure_subplot_type:
    """Wrapper that creates a simple plot using subplot routines and saves it.

    Args:
        fpath (str): The filepath for the saved plot.
        subplot_func (func): The subplot function.
        args (list): List of arguments for the subplot funciton.
        figsize (optional, tuple): The size of the figure to make.

    Returns:
        fig (matplotlib.figure.Figure): Figure for plotting.
        ax (matplotlib.axes._subplots.AxesSubplot): Axes for plotting.
    """
    # Prep the figure.
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Plot the subplot
    subplot_func(ax, *args)

    # Save and return the figure.
    if fpath is not None:
        fig.savefig(fpath, dpi=400, facecolor="w", bbox_inches="tight")

    return fig, ax
