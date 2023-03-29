import pytest
import numpy as np
import pandas as pd

import disruptivity


# PyTest Fixture that makes the data frames and index dicts.
@pytest.fixture
def gen_dataframe_and_indices():
    # Create the data frame
    CONFIG = {"elem1": [0] + [1] * 5 + [2] * 3 + [3]}
    dataframe = pd.DataFrame.from_dict(CONFIG)

    # Indices - 1x0, 2x1, 2x2, 1x3
    indicies = [0, 1, 5, 6, 7, 9]

    return dataframe, indicies


def test_indices_to_histogram1d(gen_dataframe_and_indices):
    # Unpack the fixture
    dataframe, indices = gen_dataframe_and_indices

    # Compute the histogram
    x_range = [-0.5, 3.5]

    # Compute the histogram
    hist = disruptivity.indices_to_histogram1d(
        dataframe["elem1"], x_range, indices, x_nbins=4
    )

    # Asserts
    assert (hist[0] == [1, 2, 2, 1]).all(), "Incorrect Histogram Calculation"
    assert (
        hist[1] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"


def test_indices_to_histogram2d(gen_dataframe_and_indices):
    # Unpack the fixture
    dataframe, indices = gen_dataframe_and_indices

    # Compute the histogram
    x_range = [-0.5, 3.5]

    # Compute the histogram
    hist = disruptivity.indices_to_histogram2d(
        dataframe["elem1"],
        dataframe["elem1"],
        x_range,
        x_range,
        indices,
        x_nbins=4,
        y_nbins=4,
    )

    result = np.diag([1, 2, 2, 1])

    # # Asserts
    assert (hist[0] == result).all(), "Incorrect Histogram Calculation"
    assert (
        hist[1] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"
    assert (
        hist[2] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"


def test_compute_disruptivity():
    # Bins
    num_val = np.array([0, 1, 2, 1, 0])
    denom_val = np.array([1, 1, 0, 1, 1])
    bins = np.array([0, 1, 2, 3, 4, 5])
    disrupt_result = np.array([0, 1, 0, 1, 0])
    error_result = np.array([0, np.sqrt(2), 0, np.sqrt(2), 0])
    dt = 1

    # Test the 1D Case
    # Numerator and Denominator 1D
    num1d = [num_val, bins]
    denom1d = [denom_val, bins]
    result1d = [disrupt_result, error_result]

    # Compute
    disrupt1d = disruptivity.compute_disruptivity(num1d, denom1d, dt)

    # Test
    assert (
        result1d[0] == disrupt1d[0]
    ).all(), "Error in 1D Disruptivity: Bad Disruptivity."
    assert (
        result1d[1] == disrupt1d[1]
    ).all(), "Error in 1D Disruptivity: Bad Error."

    # Test the 2D Case
    # Numerator and Denominator 2D
    num2d = [np.diag(num_val), bins, bins]
    denom2d = [np.diag(denom_val), bins, bins]
    result2d = [np.diag(disrupt_result), np.diag(error_result)]

    # Compute
    disrupt2d = disruptivity.compute_disruptivity(num2d, denom2d, dt)

    # Test
    assert (
        result2d[0] == disrupt2d[0]
    ).all(), "Error in 2D Disruptivity: Bad Disruptivity."
    assert (
        result2d[1] == disrupt2d[1]
    ).all(), "Error in 2D Disruptivity: Bad Error."
