"""probability.py contains computation routines for converting disruptivity calculations to
probabilities of disruptions."""

# imports
import numpy as np


def p_data(d: float | np.ndarray, dt: float | np.ndarray, is_disrupt: bool):
    """Computes the probability of surviving/disrupting in the next time interval dt given d.

    Args:
        d (float|np.ndarray): The disruptivity or array of disruptivities.
        dt (float|np.ndarray): The time interval for the Poisson statistics.
        is_disrupt (bool): Compute the probability for disruption or survival.

    Returns:
        [float|np.ndarray]: The probability of the the data disrupting or surviving.
    """

    if is_disrupt:
        return 1 - np.exp(-dt * d)
    return np.exp(-dt * d)


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

    # Safety
    assert d >= 0, "Invalid disruptivity, must be greater than 0."

    # Compute the negative log likelihood
    neg_log_likelihood = -np.sum(np.log(1 - np.exp(-dis_dt * d))) - np.sum(
        (-non_dis_dt * d)
    )

    return neg_log_likelihood
