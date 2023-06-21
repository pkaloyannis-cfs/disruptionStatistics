"""# Description

This file contains data loading routines to create the pandas dataframe
and index dictionaries.

Note: Nothing in this file is tested.

# Functions
"""

# Imports
import numpy as np
from scipy.io import loadmat
import pandas as pd


def load_disruptions_mat(filepath: str):
    """Load an SQL table as a .mat file. This will most likely be deprecated as new loading
    functions become available and as the format for different machines get added than CMOD.

    Args:
        filepath (str): The filepath to the .mat file to be loaded.

    Returns:
        data_df (type): Pandas dataframe of the non index data.
        index_dict (dict): Dictionary of lists of indices for indexing dataframe.
    """

    # Load the mat file
    data = loadmat(filepath)

    # Remove .mat file metadata
    data.pop("__header__")
    data.pop("__version__")
    data.pop("__globals__")

    # Skipped Keys -- For DIII-D Data
    skipped_keys = ["iperr", "ipprog"]

    # Get the Key list
    key_list = list(data.keys())

    # Flatten the arrays so that pandas can create a dataframe from it
    # and extract the index keys into their own dictionary.
    # Note: this is not the most efficient code, but it is compact
    index_dict = {}
    for key in key_list:
        # Data flattening for Pandas.
        # Offset by 1 to convert from
        # MATLAB indexing :vomit:
        data[key] = data[key].flatten()

        # Index Extraction
        if "indices" in key:
            index_dict[key] = data.pop(key) - 1

        # Skip the skipped keys
        if key in skipped_keys:
            data.pop(key)

    # Create the Pandas DataFrame
    data_df = pd.DataFrame.from_dict(data)

    # Print Load Information
    n_shots = np.unique(data_df.shot).shape[0]
    n_shots_no_disrupt = np.unique(
        data_df.shot[index_dict["indices_no_disrupt"]]
    ).shape[0]
    n_shots_disrupt = np.unique(
        data_df.shot[index_dict["indices_disrupt"]]
    ).shape[0]
    assert n_shots_disrupt + n_shots_no_disrupt == n_shots, (
        "Number of disrupts plus number of non disruptions does not equal the"
        " total shot number"
    )
    print(
        f"Total Shot Number: {n_shots}, Non-Disrupted Shots:"
        f" {n_shots_no_disrupt}, Disrupted Shots: {n_shots_disrupt}"
    )

    # Return
    return (data_df, index_dict)
