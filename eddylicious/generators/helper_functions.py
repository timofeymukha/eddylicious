import numpy as np
from scipy.integrate import simps

__all__ = ["blending_function", "delta_99", "delta_star", "theta",
           "chunks_and_offsets", "blending_function_theta"]


# Blending function for inner and outer scales as defined by Lund et al
def blending_function(eta, alpha=4, b=0.2):
    """Return the value of the blending function W for Lund's rescaling.

    Return the value of the blending function for inner and
    and outer profiles produced by Lund's rescaling.
    For eta>1 the function returns 1.

    Parameters
    ----------
    eta : ndarray
        The values of the non-dimensionalized wall-normal
        coordinate.
    alpha : float, optional
        The value of alpha (default 4.0).
    b : float
        The value of b (default 0.2).

    """
    val = np.zeros(eta.shape)

    for i in range(eta.size):
        if eta[i] <= 1:
            val[i] = 0.5*(1+1/np.tanh(alpha)*np.tanh(alpha*(eta[i]-b) /
                          ((1-2*b)*eta[i]+b)))
        else:
            val[i] = 1.0
    return val


def blending_function_theta(eta, alpha=15, b=2):
    """Return the value of the blending function for matching inner and
    outer profiles.

    Return the value of the blending function for inner and and outer
    profiles produced by Lund's rescaling. Suitable for values of eta
    obtained when momentum thickness is used as the outer scale. For
    eta>10 the function returns 1.

    Parameters
    ----------
    eta : ndarray
        The values of the non-dimensionalized wall-normal
        coordinate.
    alpha : float, optional
        The value of alpha (default 15.0).
    b : float
        The value of b (default 2).

    """
    val = np.zeros(eta.shape)

    for i in range(eta.size):
        if eta[i] <= 10:
            val[i] = 0.5*(1+1/np.tanh(alpha)*np.tanh(alpha*(eta[i]-b) /
                          ((10-2*b)*eta[i]+b)))
        else:
            val[i] = 1.0
    return val


def delta_99(y, v):
    """Compute :math:`\delta_{99}`."""

    for i in range(y.size):
        if v[i] >= 0.99*np.max(v):
            delta99 = y[i-1]
            break
    return delta99


def delta_star(y, v):
    """Compute the displacement thickness :math:`\delta^*` using
    Simpson's method."""

    return simps((1-v/v[-1]), x=y)


def theta(y, v):
    """Compute the momentum thickness using Simpson's method."""

    return simps(v/v[-1]*(1-v/v[-1]), x=y)


def chunks_and_offsets(nProcs, size):
    """Given the size of a 1d array and the number of
    processors, compute chunk-sizes for each processor
    and the starting indices (offsets) for each 
    processor.

    Parameters
    ----------
    nProcs : int
        The amount of processors.
    size : int
        The size of the 1d array to be distributed.

    Returns
    -------
        List of two 1d ndarrays of size nProcs.
        The first array contains the chunk-size for each
        processor.
        The second array contains the offset (starting 
        index) for each processor.
    """

    # To ensure integer division later
    nProcs = int(nProcs)

    try:
        # Number of procs should be positive
        assert nProcs > 0

        # All procs should have at least a 1-sized chunk
        assert nProcs <= size
    except AssertionError as e:
        raise AssertionError(e.message +
                             "Number of processors (", nProcs, ") is invalid.")

    chunks = np.zeros(nProcs, dtype=np.int64)
    nrAlloced = 0

    for i in range(nProcs):
        remainder = size - nrAlloced
        buckets = (nProcs - i)
        chunks[i] = remainder / buckets
        nrAlloced += chunks[i]

    # Calculate the offset for each processor
    offsets = np.zeros(chunks.shape, dtype=np.int64)

    for i in range(offsets.shape[0]-1):
        offsets[i+1] = np.sum(chunks[:i+1])

    try:
        assert np.sum(chunks) == size
    except AssertionError as e:
        raise AssertionError(e.message + "Chunks don't sum up to array size.")


    return [chunks, offsets]
