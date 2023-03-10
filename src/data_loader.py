"""dataloader.py is houses dataloading functions."""

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

    # Create the Pandas DataFrame
    data_df = pd.DataFrame.from_dict(data)

    # Return
    return (data_df, index_dict)
