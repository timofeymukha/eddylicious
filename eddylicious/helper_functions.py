# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

""" Various helper functions used by the generators."""

import numpy as np
from scipy.integrate import simps
from mpi4py import MPI
import os
import h5py

__all__ = ["blending_function", "delta_99", "delta_star", "momentum_thickness",
           "chunks_and_offsets", "config_to_dict", "set_write_path"]


def blending_function(eta, alpha=4, b=0.2):
    """Return the value of the blending function W for Lund's rescaling.

    Return the value of the blending function for the inner and
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


def delta_99(y, v):
    """Compute :math:`\delta_{99}`.

    Parameters
    ----------
    y : ndarray
        The independent variable.
    v : ndarray
        The velocity values.

    Returns
    -------
    float
        The value of :math:`\delta_{99}`.

    """
    delta99 = 0.0
    for i in range(y.size):
        if v[i] > 0.99*v[-1]:
            delta99 = y[i-1]
            break

    if delta99 == 0:
        raise ValueError("Error while computing delta99.")

    return delta99


def delta_star(y, v):
    """Compute the displacement thickness using Simpson's method.

    Parameters
    ----------
    y : ndarray
        The independent variable.
    v : ndarray
        The velocity values.

    Returns
    -------
    float
        The value of the displacement thickness.

    """
    return simps((1-v/v[-1]), x=y)


def momentum_thickness(y, v):
    """Compute the momentum thickness using Simpson's method.

    Parameters
    ----------
    y : ndarray
        The independent variable.
    v : ndarray
        The velocity values.

    Returns
    -------
    float
        The value of the momentum thickness.

"""
    return simps(v/v[-1]*(1-v/v[-1]), x=y)


def chunks_and_offsets(nProcs, size):
    """Given the size of a 1d array and the number of processors,
    compute chunk-sizes for each processor and the starting indices
    (offsets) for each processor.

    Parameters
    ----------
    nProcs : int
        The amount of processors.
    size : int
        The size of the 1d array to be distributed.

    Returns
    -------
    List of two ndarrays.
        The first array contains the chunk-size for each processor.
        The second array contains the offset (starting index) for
        each processor.
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


def config_to_dict(configFile):
    """Parse a config file to a dictionary."""

    configDict = {}

    for line in configFile:
        if (line[0] == '#') or (line == '\n'):
            continue
        configDict[line.split()[0]] = line.split()[1]

    return configDict


def set_write_path(config):
    """Set the writePath variable in concordance with the writer.

    For the ofnative writer: the path to constant/boundaryData directory.
    For the hdf5 writer: the hdf5 file itself.

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    writer = config["writer"]
    writePath = config["writePath"]

    if writer == "ofnative":
        inletPatchName = config["inletPatchName"]
        writePath = os.path.join(writePath, "constant", "boundaryData",
                                 inletPatchName)
        if rank == 0:
            if not os.path.exists(writePath):
                os.makedirs(writePath)

    elif writer == "hdf5":
        writePath = os.path.join(writePath, config["hdf5FileName"])
        # If the hdf5 file exists, delete it.
        if rank == 0 and os.path.isfile(writePath):
            print("HDF5 database already exists. It it will be overwritten.")
            os.remove(writePath)

        # We change the writePath to be the hdf5 file itself
        writePath = h5py.File(writePath, 'a', driver='mpio',
                              comm=MPI.COMM_WORLD)
    elif writer == "vtk":
        pass
    else:
        raise ValueError("Unknown writer: "+writer)

    return writePath
