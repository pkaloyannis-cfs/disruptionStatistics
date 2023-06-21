"""disruptivity.py contains computation routines for computing disruptivity maps."""

import pandas as pd
import numpy as np
import scipy
from scipy.stats import binned_statistic_dd
import warnings

hist1D_type = tuple[np.ndarray, np.ndarray]
hist2D_type = tuple[np.ndarray, np.ndarray, np.ndarray]

# TODO: Entry dictionary docs
# TODO: Fix the return docs on disruptivity


def indices_to_histogram(
    dataframe: pd.core.frame.DataFrame,
    entry_dict: dict,
    indices: np.ndarray,
    nbins=25,
) -> scipy.stats._binned_statistic.BinnedStatisticddResult:
    """Generate histogram from dataframe and indices over a general number of axes.

    Args:
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        entry_dict (dict): The entry dictionary (docs coming soonTM).
        indices (np.ndarray): The list of indices of note.
        nbins (int, optional): The number of histogram bins. Non square bins
            are untested so only ints are supported for now. Defaults to 25.

    Returns:
        scipy.stats._binned_statistic.BinnedStatisticddResult: The output histogram.
    """

    # Make use of dictionary ordering as of Python 3.7
    entry_range_list = []
    entry_data_list = []
    for entry in entry_dict:
        # First, select a mask value that is out of the range
        entry_range = entry_dict[entry]["range"]
        entry_range_list.append(entry_range)

        # Next, get the data from the dataframe and convert it to
        # an array. Mask the non-finite values with the mask value.
        # The mask is defined as 1.1x the maximum range parameter
        entry_data = dataframe[entry][indices]
        entry_data[~np.isfinite(dataframe[entry])] = 1.1 * entry_range[1]
        entry_data_list.append(entry_data)

    hist_dd = binned_statistic_dd(
        entry_data_list,
        None,
        range=entry_range_list,
        bins=nbins,
        expand_binnumbers=True,
        statistic="count",
    )

    return hist_dd


def compute_variable_time(
    dataframe: pd.core.frame.DataFrame,
    denom_dd: scipy.stats._binned_statistic.BinnedStatisticddResult,
    denom_indices: np.ndarray,
    nbins=25,
) -> np.ndarray:
    """Compute the variable time array for variable time disruptivity.

    Args:
        dataframe (pd.core.frame.DataFrame): The tokamak dataframe
        denom_dd (scipy.stats._binned_statistic.BinnedStatisticddResult): The denominator histogram.
        denom_indices (np.ndarray): The denominator indices.
        nbins (int, optional): Number of histogram bins. Defaults to 25.

    Returns:
        np.ndarray: The net time spent in each histogram bin.
    """

    # Dimension Calculation
    dimension = len(denom_dd.bin_edges)
    if dimension == 1:
        binnumber = np.array([denom_dd.binnumber])
    else:
        binnumber = denom_dd.binnumber
    # The : indexes all dims, : indexes the data index, np.newaxis
    # prepares the data for np.ix_() to be called
    binnumber = binnumber[:, :, np.newaxis]

    dt_array = np.zeros(denom_dd.statistic.shape)
    # test_array = np.zeros(n_disrupt.shape)
    for i in range(denom_dd.binnumber.shape[-1]):
        # Get the entry and make sure it is not out of bounds
        # For details on this, refer to the notes section of the
        # stats.binned_statistic_dd function. Basically bin 0
        # and bin nbin+1 are padded bins for outside boundaries
        entry = binnumber[:, i] - 1

        # Check if the data is out of bounds
        if (entry < 0).any() or (entry > nbins - 1).any():
            continue

        # Convert entry into a struct for array indexing.
        entry = np.ix_(*entry)

        # Get the dataframe entry for this histogram entry
        dataframe_index = denom_indices[i]
        # test_array[entry]+=1

        # Check out of bounds
        if dataframe_index - 1 < 0:
            continue
        # Check if last frame in the shot
        elif (
            dataframe["shot"][dataframe_index]
            != dataframe["shot"][dataframe_index - 1]
        ):
            continue
        # If all passes, compute dt and add it.
        else:
            dt_array[entry] += (
                dataframe["time"][dataframe_index]
                - dataframe["time"][dataframe_index - 1]
            )

    return dt_array


def compute_disruptivity(
    dataframe: pd.core.frame.DataFrame,
    entry_dict: dict,
    num_indices: np.ndarray,
    denom_indices: np.ndarray,
    nbins=25,
    dt=None,
):
    """Compute the disruptivity for certain indices of shots.

    Args:
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        entry_dict (dict): The entry dictionary (docs coming soonTM).
        num_indices (np.ndarray): The list of numerator (disruption) indices of note.
        denom_indices (np.ndarray): The list of denominator (total) indices of note.
        nbins (int, optional): The number of histogram bins. Non square bins
            are untested so only ints are supported for now. Defaults to 25.
        dt (float, optional): A fixed timestep for fixed timestep calculations.
    """
    # Compute the histograms
    num_dd = indices_to_histogram(dataframe, entry_dict, num_indices, nbins)
    denom_dd = indices_to_histogram(dataframe, entry_dict, denom_indices, nbins)

    # Parse the data
    n_disrupt = num_dd.statistic
    n_total = denom_dd.statistic

    # Compute dt using the denominator histogram
    if dt is None:
        dt_array = compute_variable_time(
            dataframe, denom_dd, denom_indices, nbins
        )
    else:
        # Fixed Timestep Computation
        dt_array = dt * n_total

    # Surpress printing Division by 0 warnings since we handle them manually.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Disruptivity equation sourced from verbal description in deVries 2009.
        # Here we handle divisions by 0 by setting the division result to 0.
        # Note that this is not the exact equation from the paper. Since
        # dt_array is a sum of all the times, not an average, we can just
        # divide by it instead.
        disruptivity = n_disrupt / dt_array
        disruptivity[~np.isfinite(disruptivity)] = 0

        # Error calculation assumes histogram error e = sqrt(n)
        # and propagates it over a division.
        # Here we handle divisions by 0 by setting the division result to 0.
        error = disruptivity * np.sqrt(1 / n_disrupt + 1 / n_total)
        error[~np.isfinite(error)] = 0

    return disruptivity, error, num_dd.bin_edges, entry_dict


# ---- Indexing ----
# Get the possible warning times
def get_indices_disruptivity(
    tokamak_config: dict,
    dataframe: pd.core.frame.DataFrame,
    index_dict: dict,
    tau=0,
    window=1,
    shotlist=None,
) -> dict:
    """Adds the indices for detectable disruptivity to the index dictionary.

    Args:
        tokamak_config (dict): The Tokamak Config Dictionary.
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        index_dict (dict): The Index Dictionary.
        tau (float)(optional): The tau class of the indices in ms. Default is 0 ms.
        window (float)(optional): The window size for data selection on either side of the data in ms (tau-window, tau). Default is 1ms.
        shotlist (list)(optional): The list of shot numbers to check for disruptions in. Default is None for no filtering.

    Returns:
        dict: The modified index dictionary.
    """

    # Find Viable Times for Detectable Disruptivity and Denominator.
    time_until_disrupt_ms = dataframe.time_until_disrupt * 1000
    tau_offset_times = time_until_disrupt_ms - tau
    warning_times = np.logical_and(
        tau_offset_times <= window, tau_offset_times >= 0
    )

    # TODO SHOT NUMBER FILTERING
    disruptive_indices = index_dict["indices_flattop_disrupt_in_flattop"]
    if shotlist is not None:
        shotlist_bool = np.isin(dataframe.shot, shotlist)
        disruptive_indices = dataframe[shotlist_bool].index

    # Assert if there is a non-disrupted pulse in the shot list
    assert (
        True
    ), f"Pulse(s) {None} in the shot list overlap with non-disruptive indices."

    # Make sure they are part of our disruptions of interest
    indices_n_disrupt = np.array(
        np.intersect1d(
            disruptive_indices,
            dataframe.loc[warning_times].index,
        )
    )

    # Create the denominator
    indices_n_total = np.append(
        index_dict["indices_flattop_no_disrupt"], disruptive_indices
    )

    return indices_n_disrupt, indices_n_total
