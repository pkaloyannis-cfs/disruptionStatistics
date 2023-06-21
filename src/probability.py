"""probability.py contains computation routines for converting disruptivity calculations to
probabilities of disruptions."""

# imports
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize, LinearConstraint


# Create a modified variable dt function that stores the sequence of dts and marks them as disrupted or not.
def compute_dt_bin(
    dataframe: pd.core.frame.DataFrame,
    denom_dd: scipy.stats._binned_statistic.BinnedStatisticddResult,
    denom_indices: np.ndarray,
    tau=0,
    window=1,
    shotlist=None,
):
    """Computes the sequence lengths dt for each histogram bin.

    Args:
        dataframe (pd.core.frame.DataFrame): The tokamak dataframe
        denom_dd (scipy.stats._binned_statistic.BinnedStatisticddResult): The denominator histogram.
        denom_indices (np.ndarray): The denominator indices.
        nbins (int, optional): Number of histogram bins. Defaults to 25.
        tau (TODO)
        window (TODO)

    Returns:
        TODO TODO TODO
    """

    # Dimension Calculation
    dim = len(denom_dd.bin_edges)
    n_bins = len(denom_dd.bin_edges[0]) - 1
    if dim == 1:
        binnumber = np.array([denom_dd.binnumber])
    else:
        binnumber = denom_dd.binnumber
    # The : indexes all dims, : indexes the data index, np.newaxis
    # prepares the data for np.ix_() to be called
    binnumber = binnumber[:, :, np.newaxis]

    dt_array = np.zeros(denom_dd.statistic.shape)

    # Store the last entry for continuity checks.
    # Here we assume that subsequent dataframe indices
    # that are in the same bin must belong to the same
    # pulse.
    #
    # Thus we must check that the last entry and last
    # dataframe index number are one and the same.
    last_entry = binnumber[:, 0] - 1
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

    # Counter for what memory location to use
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

    return pd.DataFrame.from_dict(pd_dict), max_counter


def entry_list_to_arrays(entry_df, max_counter, hist):
    """Convert the entry list of dt slices into arrays for computations.

    Args:
        entry_df (_type_): The entry dataframe produced by compute_dt_bin.
        max_counter (_type_): _description_
        hist (_type_): _description_

    Returns:
        _type_: _description_
    """

    # We do this in 1D
    non_disrupt = entry_df[entry_df.is_disrupt == False]
    disrupt = entry_df[entry_df.is_disrupt == True]

    # Preallocate the array we need
    bin_edges = np.array(hist.bin_edges)
    dim = bin_edges.shape[0]
    n_bins = bin_edges.shape[1] - 1

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
    dis_dt: np.ndarray, non_dis_dt: np.ndarray, guess=1, min_d=1e-4
) -> np.ndarray:
    """_summary_

    Args:
        dis_dt (np.ndarray): _description_
        non_dis_dt (np.ndarray): _description_
        guess (int, optional): _description_. Defaults to 1.

    Returns:
        np.ndarray: _description_
    """

    # Create the array that will store the disruptivity values
    d_array = np.zeros(dis_dt.shape[:-1])

    # Create the optimization constraints
    constraint = LinearConstraint([1], lb=min_d)

    # nditer over that value, use the multi index to index the dts
    with np.nditer(
        d_array, flags=["multi_index"], op_flags=["readwrite"]
    ) as it:
        for d in it:
            # Iterator index
            index = it.multi_index

            # If the bin is empty, continue without optimizing
            if dis_dt[index][0] == 0 and non_dis_dt[index][0] == 0:
                continue

            # Minimization routine
            res = minimize(
                disruptivity_neg_log_likelihood,
                x0=guess,
                args=(dis_dt[index], non_dis_dt[index]),
                constraints=(constraint),
            )

            # Don't draw the contraints
            if res.x[0] <= 2 * min_d:
                continue

            # Set the disruptivity
            d[...] = res.x[0]
    return d_array
