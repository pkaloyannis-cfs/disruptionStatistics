"""disruptivity.py contains computation routines for computing disruptivity maps."""

import pandas as pd
import numpy as np
import warnings

hist1D_type = tuple[np.ndarray, np.ndarray]
hist2D_type = tuple[np.ndarray, np.ndarray, np.ndarray]


def indices_to_histogram1d(
    x_series: pd.core.series.Series,
    x_range: list,
    indices: np.ndarray,
    x_nbins=25,
) -> hist1D_type:
    """Converts a pandas series and indices to a numpy 1D histogram.

    Args:
        x_series (pd.core.series.Series): The pandas series of data to be converted.
        x_range (list): The parameter limits for the histogram [min,max].
        indices (np.ndarray): The array of indices to index the pandas series.
        x_nbins (int, optional): Number of histogram bins. Defaults to 25.

    Returns:
        hist (hist1D_type): Numpy 1D histogram (histogram, bins).
    """

    # Create histrogram bins
    bins = np.histogram_bin_edges(x_series, bins=x_nbins, range=x_range)
    hist = np.histogram(x_series[indices], bins)

    return hist


def indices_to_histogram2d(
    x_series: pd.core.series.Series,
    y_series: pd.core.series.Series,
    x_range: list,
    y_range: list,
    indices: np.ndarray,
    x_nbins=25,
    y_nbins=25,
) -> hist2D_type:
    """Converts two pandas series and indices to a numpy 2D histogram.

    Args:
        x_series (pd.core.series.Series): The pandas x series of data to be converted.
        y_series (pd.core.series.Series): The pandas y series of data to be converted.
        x_range (list): The parameter x limits for the histogram [min,max].
        y_range (list): The parameter y limits for the histogram [min,max].
        indices (np.ndarray): The array of indices to index the pandas series.
        x_nbins (int, optional): Number of x histogram bins. Defaults to 25.
        y_nbins (int, optional): Number of y histogram bins. Defaults to 25.

    Returns:
        hist (hist2D_type): Numpy 2D histogram (histogram, x_bins, y_bins).
    """

    hist = np.histogram2d(
        x_series[indices],
        y_series[indices],
        bins=[x_nbins, y_nbins],
        range=[x_range, y_range],
    )
    return hist


def compute_disruptivity(
    hist_disrupt,
    hist_total,
    dt: float,
):
    """Computes the disruptivity and errorbars.

    Args:
        hist_disrupt : Numerator of disruption events.
        hist_total : Denominator of all data points.
        dt (float): Time step (s).

    Returns:
        disruptivity : The disruptivity.
        error : The errorbars.
    """

    # Parse the data
    n_disrupt = hist_disrupt[0]
    n_total = hist_total[0]

    # Asserts
    assert len(n_disrupt) == len(n_total), (
        "Numerator and Denominator array size mismatch. Numerator size:"
        f" {len(n_disrupt)}, Denominator size: {len(n_total)}"
    )
    assert dt != 0, "dt is 0."
    assert len(hist_disrupt) == len(hist_total), (
        "Numerator and Denominator have different dimensions."
        f"Numerator Dimension {len(hist_disrupt)-1}, "
        f"Denominator Dimension {len(hist_disrupt)-1}, "
    )

    # Surpress printing Division by 0 warnings since we handle them manually.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Disruptivity equation sourced from verbal description in deVries 2009.
        # Here we handle divisions by 0 by setting the division result to 0.
        disruptivity = n_disrupt / n_total / dt
        disruptivity[~np.isfinite(disruptivity)] = 0

        # Error calculation assumes histogram error e = sqrt(n)
        # and propagates it over a division.
        # Here we handle divisions by 0 by setting the division result to 0.
        error = disruptivity * np.sqrt(1 / n_disrupt + 1 / n_total)
        error[~np.isfinite(error)] = 0

    return (disruptivity, error)
