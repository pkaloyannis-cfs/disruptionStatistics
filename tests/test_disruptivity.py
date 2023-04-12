import pytest
import numpy as np
import pandas as pd

import disruptivity


# PyTest Fixture that makes the data frames and index dicts.
@pytest.fixture
def gen_dataframe_and_indices():
    # Create the data frame
    CONFIG = {
        "elem1": [0] + [1] * 5 + [2] * 3 + [3],
        "elem2": [0] + [1] * 5 + [2] * 3 + [3],
        "time": np.arange(10),
        "shot": [1] * 10,
    }

    dataframe = pd.DataFrame.from_dict(CONFIG)

    # Indices - 1x0, 2x1, 2x2, 1x3
    indices = [0, 1, 5, 6, 7, 9]

    # Create the entry_dict
    entry_dict = {
        "elem1": {"range": [-0.5, 3.5]},
        "elem2": {"range": [-0.5, 3.5]},
    }

    return dataframe, indices, entry_dict


def test_indices_to_histogram(gen_dataframe_and_indices):
    # Unpack the fixture
    dataframe, indices, entry_dict = gen_dataframe_and_indices

    # Compute the histogram general dimension case
    hist_2d = disruptivity.indices_to_histogram(
        dataframe, entry_dict, indices, nbins=4
    )

    result = np.diag([1, 2, 2, 1])

    # # Asserts
    assert (
        hist_2d.statistic == result
    ).all(), "Incorrect Histogram Calculation"
    assert (
        hist_2d.bin_edges[0] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"
    assert (
        hist_2d.bin_edges[1] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"

    # The special 1D case.
    entry_dict.pop("elem2")
    # Compute the histogram general dimension case
    hist_1d = disruptivity.indices_to_histogram(
        dataframe, entry_dict, indices, nbins=4
    )

    # Asserts
    assert (
        hist_1d.statistic == [1, 2, 2, 1]
    ).all(), "Incorrect Histogram Calculation"
    assert (
        hist_1d.bin_edges[0] == [-0.5, 0.5, 1.5, 2.5, 3.5]
    ).all(), "Incorrect Bin Calculation"


def test_compute_variable_time_array(gen_dataframe_and_indices):
    # Unpack the fixture and prep variables
    dataframe, indices, entry_dict = gen_dataframe_and_indices
    nbins = 4
    denom_dd = disruptivity.indices_to_histogram(
        dataframe, entry_dict, indices, nbins
    )

    # The first test is the mapping test, that is
    # given a denom histogram, if dt=1 everywhere
    # match the statistics.
    dt_array = disruptivity.compute_variable_time(
        dataframe, denom_dd, indices, nbins
    )

    # Fill the first time slice, since no dt can be computed for it
    dt_array[0, 0] += 1

    assert (
        dt_array == denom_dd.statistic
    ).all(), "Incorect Variable Timestep Mapping or Calculation"

    # The next test is the out of bounds test
    # When histogram entries are out of bounds
    # They should be ignored.
    dataframe["elem1"][0] = -1
    denom_dd = disruptivity.indices_to_histogram(
        dataframe, entry_dict, indices, nbins
    )
    dt_array = disruptivity.compute_variable_time(
        dataframe, denom_dd, indices, nbins
    )
    dataframe["elem1"][0] = 0

    # Don't need to fill the first time slice since it is out of bounds.

    assert (dt_array == denom_dd.statistic).all(), "Incorrect Boundary Handling"

    # The next test is shot number transition test.
    # On a shot number transition, we expect to ignore
    # that dt since it cannot be computed.
    dataframe["shot"] = [1] * 6 + [2] * 4
    denom_dd = disruptivity.indices_to_histogram(
        dataframe, entry_dict, indices, nbins
    )
    dt_array = disruptivity.compute_variable_time(
        dataframe, denom_dd, indices, nbins
    )

    # Fill the first time slice, since no dt can be computed for it
    dt_array[0, 0] += 1

    # Fill the [2,2] index corresponding to missing 2 in the summation.
    # that comes from the pulse number transition
    dt_array[2, 2] += 1

    assert (
        dt_array == denom_dd.statistic
    ).all(), "Incorrect Shot Transition Handling"

    return


def test_compute_disruptivity(gen_dataframe_and_indices):
    # Unpack the fixture
    dataframe, indices, entry_dict = gen_dataframe_and_indices

    # 1x0 5x1 3x2 1x3
    indices_2 = np.arange(0, 10)

    # The answers to the fixed time tests
    num_val = np.array([1, 2, 2, 1])
    denom_val = np.array([1, 5, 3, 1])
    disrupt_result = num_val / denom_val
    error_result = disrupt_result * np.sqrt(
        [2, 1 / 2 + 1 / 5, 1 / 2 + 1 / 3, 2]
    )
    dt = 1

    entry_dict.pop("elem2")

    # Compute the fixed time 1D dirsuptivity
    disrupt1d, errors1d, bin_edges_1d, entry_dict = (
        disruptivity.compute_disruptivity(
            dataframe, entry_dict, indices, indices_2, nbins=4, dt=dt
        )
    )

    # Test
    assert (
        disrupt1d == disrupt_result
    ).all(), "Error in 1D Disruptivity: Bad Disruptivity."
    assert (
        errors1d == error_result
    ).all(), "Error in 1D Disruptivity: Bad Error."
