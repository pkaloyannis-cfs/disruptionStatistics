"""disruptivity.py contains computation routines for computing disruptivity maps."""

import pandas as pd
import numpy as np
import scipy
from scipy.stats import binned_statistic_dd
from scipy.optimize import minimize, LinearConstraint
import warnings

hist1D_type = tuple[np.ndarray, np.ndarray]
hist2D_type = tuple[np.ndarray, np.ndarray, np.ndarray]

# TODO: Entry dictionary docs

# Data masks for different limiting cases
# This allows special processing of these
# Data points later for interpolation or plotting.
NODATAMASK = -1
NODISRUPTIONMASK = -2
ALLDISRUPTIONMASK = -3


# ---- Sampling Method Functions ----
def indices_to_histogram(
    dataframe: pd.core.frame.DataFrame,
    entry_dict: dict,
    indices: np.ndarray,
    nbins=25,
) -> scipy.stats._binned_statistic.BinnedStatisticddResult:
    """Generate histogram from dataframe and indices over a general number of axes. For use in
    sampling disruptivity method.

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
    """Compute the variable time array for variable time disruptivity. For use in sampling
    disruptivity method.

    Args:
        dataframe (pd.core.frame.DataFrame): The tokamak dataframe
        denom_dd (scipy.stats._binned_statistic.BinnedStatisticddResult): The denominator histogram.
        denom_indices (np.ndarray): The denominator indices.
        nbins (int, optional): Number of histogram bins. Defaults to 25.

    Returns:
        np.ndarray: The net time spent in each histogram bin.
    """

    # Get the bin information from the histogram
    binnumber, n_bins, dimension = get_bin_info(denom_dd)

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


def compute_disruptivity_sampling(
    dataframe: pd.core.frame.DataFrame,
    entry_dict: dict,
    num_indices: np.ndarray,
    denom_indices: np.ndarray,
    nbins=25,
    dt=None,
):
    """Compute the disruptivity for certain indices of shots using the sampling method.

    When no data is available, returns NODATAMAS (default -1).
    When no disruptions are available, returns NODISRUPTIONMASK (default -2).
    When only disruptions are available, returns ALLDISRUPTIONMASK (default -3).

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

    # Get the histogram indices to mask.
    n_no_disrupt = n_disrupt == 0
    n_all_disrupt = n_disrupt == n_total
    n_no_data = n_total == 0

    # Apply the masks
    # Order matters since no data must be applied last.
    disruptivity[n_no_disrupt] = NODISRUPTIONMASK
    disruptivity[n_all_disrupt] = ALLDISRUPTIONMASK
    disruptivity[n_no_data] = NODATAMASK

    error[n_no_disrupt] = 0
    error[n_all_disrupt] = 0
    error[n_no_data] = 0

    return disruptivity, error, num_dd.bin_edges, entry_dict


# ---- Likelihood Method Functions ----


def compute_continuous_time_bin(
    dataframe: pd.core.frame.DataFrame,
    denom_dd: scipy.stats._binned_statistic.BinnedStatisticddResult,
    denom_indices: np.ndarray,
    tau=0,
    window=1,
    shotlist=None,
):
    """Computes the sequence lengths dt for each histogram bin. For use in likelihood disruptivity
    method.

    Args:
        dataframe (pd.core.frame.DataFrame): The tokamak dataframe
        denom_dd (scipy.stats._binned_statistic.BinnedStatisticddResult): The denominator histogram.
        denom_indices (np.ndarray): The denominator indices.
        nbins (int, optional): Number of histogram bins. Defaults to 25.
        tau (float)(optional): The tau class of the indices in ms. Default is 0 ms.
        window (float)(optional): The window size for data selection on either side of the data in ms (tau-window, tau). Default is 1ms.
        shotlist (list)(optional): The list of shot numbers to check for disruptions in. Default is None for no filtering.
    Returns:
        pd_df (pd.core.frame.DataFrame): The continuous dt dictionary for maximum likelihood computations.
        max_counter (int): The maximum number of entries for a histogram bin. Used for memory allocation later on.
    """

    # Get the bin information from the histogram
    binnumber, n_bins, dim = get_bin_info(denom_dd)

    # Store the last entry for continuity checks.
    # Here we assume that subsequent dataframe indices
    # that are in the same bin must belong to the same
    # pulse.
    #
    # Thus we must check that the last entry and last
    # dataframe index number are one and the same.
    #
    # Note: last_entry is set to an array of None so it
    # can be filled with the first inbounds data later
    last_entry = np.array(None)
    dt_integrator = 0
    last_dataframe_index = denom_indices[0] - 1
    last_shot = dataframe["shot"][last_dataframe_index + 1]

    # Lists that contain the for loop outputs.
    ix_list = []
    shot_list = []
    dt_list = []
    disrupt_list = []

    # Make is disrupt global for namespacing reasons
    is_disrupt = False

    # Counter for amount of data per bin
    counter = np.zeros([n_bins] * dim, dtype=int)

    # test_array = np.zeros(n_disrupt.shape)
    for i in range(denom_dd.binnumber.shape[-1]):
        # Get the entry and make sure it is not out of bounds
        # For details on this, refer to the notes section of the
        # stats.binned_statistic_dd function. Basically bin 0
        # and bin nbin+1 are padded bins for outside boundaries
        entry = binnumber[:, i] - 1

        # Check if the data is out of bounds
        if (entry < 0).any() or (entry > n_bins - 1).any():
            continue

        # Set the initial last_entry memory to be the first entry
        # that is in bounds. Otherwise if the first data entry is
        # out of bounds this will throw an indexing error.
        if (last_entry == None).any():
            last_entry = entry

        # Get the dataframe entry for this histogram entry
        dataframe_index = denom_indices[i]
        shot = dataframe["shot"][dataframe_index]

        # Data continuity check.
        # If the current entry is not the last entry
        # place all the data into memory and increment
        # the entry. Also stop integration if we hit the
        # time of disruption.
        if (
            (entry != last_entry).all()
            or dataframe_index != last_dataframe_index + 1
            or is_disrupt
        ):
            # Convert entry into a struct for array indexing.
            ix_entry = np.ix_(*last_entry)

            # List Appends
            shot_list.append(last_shot)
            ix_list.append(ix_entry)
            dt_list.append(dt_integrator)  # +tau/1000
            disrupt_list.append(is_disrupt)

            # Increment the counter for memory allocation later
            counter[ix_entry] += 1

            # Reset the integrator and is_disrupt flags.
            last_entry = entry
            last_shot = shot
            dt_integrator = 0

        # Set the is_disrupt flag:
        is_disrupt = False
        time_until_disrupt_ms = (
            dataframe.time_until_disrupt[dataframe_index] * 1000
        )
        if (
            time_until_disrupt_ms - tau <= window
            and time_until_disrupt_ms - tau >= 0
        ):
            is_disrupt = True

        # Reset the last dataframe index
        last_dataframe_index = dataframe_index

        # Check if first frame in the shot
        # Should be impossible since we exclude ramp ups
        if shot != dataframe["shot"][dataframe_index - 1]:
            continue
        # Else: integrate the time
        else:
            dt_integrator += (
                dataframe["time"][dataframe_index]
                - dataframe["time"][dataframe_index - 1]
            )

    # After the loop, convert the lists into a pandas dataframe
    pd_dict = {
        "shot": shot_list,
        "ix": ix_list,
        "dt": dt_list,
        "is_disrupt": disrupt_list,
    }

    # Get the maximum ix entries for malloc later
    max_counter = counter.max()
    pd_df = pd.DataFrame.from_dict(pd_dict)

    return pd_df, max_counter


def entry_list_to_arrays(entry_df, max_counter, hist):
    """Convert the entry list of dt slices into arrays for computations.

    Args:
        entry_df (_type_): The entry dataframe produced by compute_dt_bin.
        max_counter (_type_): _description_
        hist (_type_): _description_

    Returns:
        _type_: _description_
    """

    # Sort the data into disrupted and non-disrupted
    non_disrupt = entry_df[entry_df.is_disrupt == False]
    disrupt = entry_df[entry_df.is_disrupt == True]

    # Preallocate the array we need
    bin_edges = np.array(hist.bin_edges)
    dim = bin_edges.shape[0]
    n_bins = bin_edges.shape[1] - 1

    #! THIS IS REPEATED CODE. WE CAN REMOVE THIS WITH A DOUBLE FOR LOOP
    #! FIX BEFORE MERGE ONTO MAIN
    # Preallocate memory to store the dts
    dis_dt_array = np.zeros([n_bins] * dim + [max_counter])
    non_dis_dt_array = np.zeros([n_bins] * dim + [max_counter])

    # Counter for what memory location to use
    dis_counter = np.zeros([n_bins] * dim, dtype=int)
    non_dis_counter = np.zeros([n_bins] * dim, dtype=int)

    for index, entry in non_disrupt.iterrows():
        # Data selection
        counter = non_dis_counter[entry["ix"]]
        non_dis_dt_array[entry["ix"] + tuple(counter)] = entry["dt"]
        non_dis_counter[entry["ix"]] += 1

    for index, entry in disrupt.iterrows():
        # Data selection
        counter = dis_counter[entry["ix"]]
        dis_dt_array[entry["ix"] + tuple(counter)] = entry["dt"]
        dis_counter[entry["ix"]] += 1

    return dis_dt_array, non_dis_dt_array


def p_data(d: float | np.ndarray, dt: float | np.ndarray, is_disrupt: bool):
    """Computes the probability of surviving/disrupting in the next time interval dt given d. If
    d*dt=0 then return 1 so as to not impact the log likelihood calculation.

    Args:
        d (float|np.ndarray): The disruptivity or array of disruptivities.
        dt (float|np.ndarray): The time interval for the Poisson statistics.
        is_disrupt (bool): Compute the probability for disruption or survival.

    Returns:
        [float|np.ndarray]: The probability of the the data disrupting or surviving.
    """

    # NEED TO COMMENT THIS
    if is_disrupt:
        return np.where(d * dt == 0, 1, 1 - np.exp(-dt * d))
    return np.where(d * dt == 0, 1, np.exp(-dt * d))


def disruptivity_neg_log_likelihood(
    d: float,
    dis_dt: np.ndarray,
    non_dis_dt: np.ndarray,
) -> float:
    """Computes the negative log likelihood for a bin of disruptivity data.

    Args:
        d (float): The disruptivity of said bin.
        dis_dt (np.ndarray): The disrupted data dt windows.
        non_dis_dt (np.ndarray): The non disrupted data dt windows.

    Returns:
        float: The negative log likelihood for the data. To be minimized.
    """

    # Compute the probabilities
    p_dis = p_data(d, dis_dt, True)
    p_non_dis = p_data(d, non_dis_dt, False)

    # Compute the negative log likelihood
    neg_log_likelihood = -np.sum(np.log(p_dis), axis=-1) - np.sum(
        np.log(p_non_dis), axis=-1
    )

    return neg_log_likelihood


def find_disruptivity(
    dis_dt: np.ndarray, non_dis_dt: np.ndarray, guess=1, min_d=1e-5
) -> np.ndarray:
    """_summary_

    Args:
        dis_dt (np.ndarray): The dis_dt array.
        non_dis_dt (np.ndarray): The non_dis_dt array.
        guess (int, optional): The initial guess for the disruptivity calculation. Defaults to 1.
        min_d (float, optional): The minimum value for d that can be optimized for.

    Returns:
        np.ndarray: _description_
    """

    # Create the array that will store the disruptivity values
    disruptivity = np.zeros(dis_dt.shape[:-1])

    # Create the optimization constraints
    constraint = LinearConstraint([1], lb=min_d)

    # nditer over that value, use the multi index to index the dts
    with np.nditer(
        disruptivity, flags=["multi_index"], op_flags=["readwrite"]
    ) as it:
        for d in it:
            # Iterator index
            index = it.multi_index

            # Masking Step
            # If the bin is empty, continue without optimizing
            if dis_dt[index][0] == 0 and non_dis_dt[index][0] == 0:
                d[...] = NODATAMASK
                continue
            # If the above check fails and there are no disruptions
            # then there must non-disrupted data
            elif dis_dt[index][0] == 0:
                d[...] = NODISRUPTIONMASK
                continue
            # If the above check fails and there are no non-disruptions
            # then there must disrupted data
            elif non_dis_dt[index][0] == 0:
                d[...] = ALLDISRUPTIONMASK
                continue

            # Minimization routine
            res = minimize(
                disruptivity_neg_log_likelihood,
                x0=guess,
                args=(dis_dt[index], non_dis_dt[index]),
                constraints=(constraint),
            )

            # Don't draw the contraints
            # This should be very rare since we skip
            # data that has no disruptions
            if res.x[0] <= 2 * min_d:
                continue

            # Set the disruptivity
            d[...] = res.x[0]
    return disruptivity


def compute_disruptivity_likelihood(
    dataframe: pd.core.frame.DataFrame,
    entry_dict: dict,
    indices: np.ndarray,
    nbins=25,
    tau=0,
    window=2,
):
    """Compute the disruptivity for certain indices of shots using the sampling method.

    When no data is available, returns NODATAMAS (default -1).
    When no disruptions are available, returns NODISRUPTIONMASK (default -2).
    When only disruptions are available, returns ALLDISRUPTIONMASK (default -3).

    Args:
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        entry_dict (dict): The entry dictionary (docs coming soonTM).
        indices (np.ndarray): The list of indices to compute disruptivity over.
        nbins (int, optional): The number of histogram bins. Non square bins
            are untested so only ints are supported for now. Defaults to 25.
        tau (float)(optional): The tau class of the indices in ms. Default is 0 ms.
        window (float)(optional): The window size for data selection on either side of the data in ms (tau-window, tau). Default is 2ms.
    """
    # Compute the histogram
    hist = indices_to_histogram(dataframe, entry_dict, indices, nbins)

    # Continuity checks
    entry_list, max_counter = compute_continuous_time_bin(
        dataframe, hist, indices, tau=tau, window=window
    )

    # Data formatting for more efficient computation
    dis_dt_list, non_dis_dt_list = entry_list_to_arrays(
        entry_list, max_counter, hist
    )

    # Minimize and apply masks as needed
    disruptivity = find_disruptivity(dis_dt_list, non_dis_dt_list, guess=1)

    # TODO ERRORBARS
    error = np.zeros(disruptivity.shape)

    return disruptivity, error, hist.bin_edges, entry_dict


# ---- Indexing & Helpers ----


def get_bin_info(hist):
    """AI is creating summary for prep_bins.

    Args:
        hist (scipy.stats._binned_statistic.BinnedStatisticddResult): The histogram to get bin data from.

    Returns:
        binnumber (np.ndarrray): Formatted bin count information for indexing.
        n_bins (int): The number of bins.
        dimension (int): The dimension of the histogram.
    """

    # Dimension Calculation
    dimension = len(hist.bin_edges)
    n_bins = len(hist.bin_edges[0]) - 1
    if dimension == 1:
        binnumber = np.array([hist.binnumber])
    else:
        binnumber = hist.binnumber
    # The : indexes all dims, : indexes the data index, np.newaxis
    # prepares the data for np.ix_() to be called
    binnumber = binnumber[:, :, np.newaxis]

    return (
        binnumber,
        n_bins,
        dimension,
    )


def get_indices_disruptivity(
    tokamak_config: dict,
    dataframe: pd.core.frame.DataFrame,
    index_dict: dict,
    tau=0,
    window=2,
    shotlist=None,
) -> dict:
    """Adds the indices for detectable disruptivity to the index dictionary.

    Args:
        tokamak_config (dict): The Tokamak Config Dictionary.
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        index_dict (dict): The Index Dictionary.
        tau (float)(optional): The tau class of the indices in ms. Default is 0 ms.
        window (float)(optional): The window size for data selection on either side of the data in ms (tau-window, tau). Default is 2ms.
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

    # Shot Number Filtering
    # Ensure the disrupted pulse
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
