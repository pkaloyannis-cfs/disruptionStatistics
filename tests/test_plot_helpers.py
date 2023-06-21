import pytest
import os
import matplotlib.pyplot as plt

import vis.plot_helpers as plot_helpers


def dummy_plotting_function(ax, arg1):
    """A dummy plotter.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): Axes for plotting.
        arg1 (int): Dummy parameter for the slope of the line to plot.
    """

    # Create a line.
    x_list = [i for i in range(100)]
    y_list = [x * arg1 for x in x_list]

    ax.plot(x_list, y_list)


def test_plot_subplot():
    """Test function for plot_subplot. Checks that it can plot a dummy plot.

    Also tests the optional file saving functionality.
    """

    # Create a plot without saving
    plot_helpers.plot_subplot(
        fpath=None,
        subplot_func=dummy_plotting_function,
        args=[1.0],
    )

    # Create a plot with saving in the run directory
    # Then check that the figure exists and delete it
    # as a cleanup step.
    path = "./test_plot.png"
    plot_helpers.plot_subplot(
        fpath=path,
        subplot_func=dummy_plotting_function,
        args=[1.0],
    )
    assert os.path.exists(path), "Function plot_subplot did not write to file."
    os.remove(path)
