"""# Description

This file contains computation routines for nonlinear boundary avoidance.

Note: No function in here is unit tested.

# Functions
"""

import numpy as np
import scipy
from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter


def iterative_fill(
    disruptivity: np.ndarray,
    bins: np.ndarray,
    d_max_mult=1,
    filter_size=11,
    max_iter=100,
    rel_tol=1e-6,
    sigma=0.5,
    truncate=3,
):
    """Iterative averaging procedure Hold real data, no data, all disrupted as fixed values
    Iteratively loop: Conv2D with averaging kernel of size n. Replace the fixed values, Check for
    convergence with last iteration via relative_tol (bin wise). Check for max iterations.

    Args:
        disruptivity (np.ndarray): The disruptivity array.
        bins (np.ndarray): The bins of the disruptivity array.
        d_max_mult (int, optional): Multiplier for the maxiumum
            fill value.. Defaults to 1.
        filter_size (int, optional): Size of the filter to use for
            averaging. Only odd numbers allowed. Defaults to 11.
        max_iter (int, optional): Max number of iterations for the
            averaging routine. Defaults to 100.
        rel_tol ([type], optional): Subsequent step convergence
            tolerance. Defaults to 1e-6.
        sigma (float, optional): std of the Gaussian filter used
            after the iterative fill. Defaults to 0.5.
        truncate (int, optional): Width of the Gaussian filter used
            after the iterative fill. Defaults to 3.

    Returns:
        (np.ndarray): The filled and smoothed disruptivity map.
    """

    # TODO: WHAT IS WRONG WITH filter_size-7

    # Create the fixed points array
    fixed_points = np.copy(disruptivity)
    max_fill = disruptivity[disruptivity > 0].max() * d_max_mult
    fixed_points[disruptivity == -1] = max_fill
    fixed_points[disruptivity == -3] = max_fill

    # Create the iter_array
    iter_array = np.zeros(fixed_points.shape)
    iter_array[disruptivity != -2] = fixed_points[disruptivity != -2]

    # Create the convolutional filter
    assert filter_size % 2 == 1, "Even Filter sizes not allowed."
    n_dims = len(bins)
    ave_filter = np.ones([filter_size] * n_dims) / (filter_size**n_dims)
    pad_size = int(np.floor(filter_size / 2))

    # Iterate
    for i in range(max_iter):
        # Step 1: Convolve
        temp = np.pad(
            iter_array, pad_size, "constant", constant_values=iter_array.max()
        )
        new_iter_array = scipy.signal.fftconvolve(
            temp, ave_filter, mode="valid"
        )

        # Step 2: Replace the fixed points
        new_iter_array[disruptivity != -2] = fixed_points[disruptivity != -2]

        # Step 3: Relative Tolerance Check
        # Check for division by 0s in early cycles
        if (new_iter_array != 0).all():
            residuals = abs(new_iter_array - iter_array) / new_iter_array
            if (residuals <= rel_tol).all():
                print(f"Data filling converged after {i} iterations.")
                break
            if i == max_iter:
                print(f"Data filling failed to converged after {i} iterations.")

        # Save the new iteration
        iter_array = new_iter_array

    # Run one last smoothing operation on the data with a small Gaussian Filter
    iter_array = gaussian_filter(iter_array, sigma=sigma, truncate=truncate)

    return iter_array


def gaussian(x, sig):
    """Gaussian distribution centered on 0.

    Used to create filters.

    Args:
        x (np.ndarray): Position to evaluate at.
        sig (float): Standard deviation.

    Returns:
        (np.ndarray): G(x, sigma)
    """
    return 1 / (2 * np.pi * np.sqrt(sig)) * np.exp(-((x / sig) ** 2) / 2)


def gaussian_x(x, sig):
    """Derivative Gaussian distribution centered on 0.

    Used to create filters.

    Args:
        x (np.ndarray): Position to evaluate at.
        sig (float): Standard deviation.

    Returns:
        (np.ndarray): G_x(x, sigma)
    """
    return -x / sig**2 * gaussian(x, sig)


def compute_grads(disruptivity, bins, sig=4, filter_size=5):
    """Compute the gradient maps for the disruptivity maps.

    Currently only implemented for 2D.

    Args:
        disruptivity (np.ndarray): The disruptivity array.
        bins (np.ndarray): The bins of the disruptivity array.
        sig (float): Standard deviation. Defaults to 4.
        filter_size (int, optional): Size of the filter to use for
            averaging. Only odd numbers allowed. Defaults to 11.

    Returns:
        (np.ndarray, np.ndarray): disruptivity_x, disruptivity_y
    """
    # Get the bin spacing
    dx = bins[0][1] - bins[0][0]
    dy = bins[1][1] - bins[1][0]

    # Create the Filter
    assert filter_size % 2 == 1, "Even Filter sizes not allowed."
    x = np.linspace(-5, 5, filter_size)
    g_filter = gaussian(x, sig)
    g_x_filter = gaussian_x(x, sig)

    # Instead of doing 2 convolve 1Ds, for now construct a 2D kernel
    # To make these real derivative will probably need some rescaling by sigma.
    max_fill = disruptivity.max()
    kernel_x = np.outer(g_filter, g_x_filter)
    grad_x = (
        convolve2d(disruptivity, kernel_x.T, mode="same", fillvalue=max_fill)
        * dx
    )
    grad_y = (
        convolve2d(disruptivity, kernel_x, mode="same", fillvalue=max_fill) * dy
    )
    return grad_x, grad_y


def create_interpers(grad_x, grad_y, bins):
    """Creates the interpolators used for boundary avoidance.

    Args:
        grad_x (np.ndarray): The disruptivity_x map.
        grad_y (np.ndarray): The disruptivity_y map.
        bins (np.ndarray): The bins of the disruptivity array.

    Returns:
        (scipy.interpolate.RegularGridInterpolator, scipy.interpolate.RegularGridInterpolator):
            The x and y gradient interpolators.
    """
    # Create two new interpolators for gx gy
    # Prep the interpolator
    bin_centers = (np.array(bins)[:, 1:] + np.array(bins)[:, :-1]) / 2
    xx = np.meshgrid(*bin_centers)
    interper_x = scipy.interpolate.RegularGridInterpolator(
        bin_centers, grad_x, method="linear", bounds_error=False, fill_value=0
    )
    interper_y = scipy.interpolate.RegularGridInterpolator(
        bin_centers, grad_y, method="linear", bounds_error=False, fill_value=0
    )
    return interper_x, interper_y
