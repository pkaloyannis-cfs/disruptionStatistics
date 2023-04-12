"""indexing.py contains several useful indexing routines for selection data for disruptivity
calculations."""

import pandas as pd
import numpy as np


# Get the possible warning times
def get_indices_detectable_disruptivity(
    tokamak_config: dict, dataframe: pd.core.frame.DataFrame, index_dict: dict
) -> dict:
    """Adds the indices for detectable disruptivity to the index dictionary.

    Args:
        tokamak_config (dict): The Tokamak Config Dictionary.
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        index_dict (dict): The Index Dictionary.

    Returns:
        dict: The modified index dictionary.
    """

    # Find Viable Times for Detectable Disruptivity and Denominator.
    disrupt_warning_time_ms = tokamak_config["disrupt_warning_time_ms"]
    disrupt_warning_window_ms = tokamak_config["disrupt_warning_window_ms"]
    time_until_disrupt_ms = dataframe.time_until_disrupt * 1000
    warning_times = (
        np.abs(time_until_disrupt_ms - disrupt_warning_time_ms)
        <= disrupt_warning_window_ms
    )
    pre_warning_times = (
        time_until_disrupt_ms
        >= disrupt_warning_time_ms - disrupt_warning_window_ms
    )

    # Intersect with the flattop to get the indices.
    # Pre warning times in flat tops only.
    indices_n_pre_warning = np.array(
        np.intersect1d(
            index_dict["indices_flattop_disrupt_in_flattop"],
            dataframe.loc[pre_warning_times].index,
        )
    )

    # Make sure they are part of our disruptions of interest
    indices_n_detectable_disrupt = np.array(
        np.intersect1d(
            index_dict["indices_flattop_disrupt_in_flattop"],
            dataframe.loc[warning_times].index,
        )
    )

    # Create the denominator
    indices_n_detectable_total = np.append(
        index_dict["indices_flattop_no_disrupt"], indices_n_pre_warning
    )

    # Add them to the index_dict.
    index_dict["indices_n_pre_warning"] = indices_n_pre_warning
    index_dict["indices_n_detectable_disrupt"] = indices_n_detectable_disrupt
    index_dict["indices_n_detectable_total"] = indices_n_detectable_total

    return index_dict


# Get the possible warning times
def get_indices_disruptivity(
    tokamak_config: dict, dataframe: pd.core.frame.DataFrame, index_dict: dict
) -> dict:
    """Adds the indices for disruptivity calculation to the index dictionary.

    Args:
        tokamak_config (dict): The Tokamak Config Dictionary.
        dataframe (pd.core.frame.DataFrame): The Tokamak Dataframe.
        index_dict (dict): The Index Dictionary.

    Returns:
        dict: The modified index dictionary.
    """

    # Find the disruption indices
    indices_n_disrupt = np.array(
        np.intersect1d(
            index_dict["indices_flattop_disrupt_in_flattop"],
            dataframe.loc[dataframe["time_until_disrupt"] == 0].index,
        )
    )

    # Create the denominator
    indices_n_total = np.append(
        index_dict["indices_flattop_no_disrupt"],
        index_dict["indices_flattop_disrupt_in_flattop"],
    )

    # Add them to the index_dict
    index_dict["indices_n_disrupt"] = indices_n_disrupt
    index_dict["indices_n_total"] = indices_n_total

    return index_dict
