# This file is part of eddylicious
# (c) 2020 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for generating  inflow velocity fields consisting
of spatially and temporally uncorrelated  fluctuations with predetermined
first- and second-order statistics.

See Appendix A in

Lund T.S., Wu X., Squires K.D. Generation of turbulent inflow
data for spatially-developing boundary layer simulations.
J. Comp. Phys. 1998; 140:233-58.

"""
from __future__ import print_function
from __future__ import division
import numpy as np
from mpi4py import MPI
from ..helper_functions import chunks_and_offsets

__all__ = ["random_fluctuations_generate"]


def random_fluctuations_generate(writerFunction,
                                 inflowPoints,
                                 uMean,
                                 reynoldsStress,
                                 dt, t0, timePrecision, tEnd):
    """Generate the inflow velocity using Gaussian random variables
    manipulated to comply to given first- and second-order statistics.

    Parameters
    ----------
    writerFunction : function
        The function to use for writing the data.
    inflowPoints: ndarray
        The coordinates of the point of the inflow patch, shape (N, 3)
    uMean: ndarray
        The mean velocity values at each point of the inflow patch,
        shape (N, 3), columns corresponding to x, y, and z coordinates.
    reynoldsStress: ndarray
        The Reynolds stress tensor at each point of the inflow patch,
        shape (N, 6), columns orderd to correspond to <uu>, <uv>, <uw>,
        <vv>, <vw>, <ww>. Algorithm will crash if the resulting matrix
        is not positive-definite.
    dt : float
        The time-step to be used in the simulation. This will be used to
        associate a time-value with the produced velocity fields.
    t0 : float
        The starting time to be used in the simulation. This will be
        used to associate a time-value with the produced velocity.
    timePrecision : int
        Number of points after the decimal to keep for the time value.
    tEnd : float
        The ending time for the simulation.
    """

    # Grab info regarding parallelization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()

    nPoints = inflowPoints.shape[0]

    # Cholesky decomposition of the Reynolds stress tensor at each point
    a = np.zeros((nPoints, 3, 3))

    # Compute the decomposition at each point
    for i in range(nPoints):
        a[i, :, :] = np.linalg.cholesky(
        [[reynoldsStress[i, 0], reynoldsStress[i, 1], reynoldsStress[i, 2]],
         [reynoldsStress[i, 1], reynoldsStress[i, 3], reynoldsStress[i, 4]],
         [reynoldsStress[i, 2], reynoldsStress[i, 4], reynoldsStress[i, 5]]])

    print(a)

    # Get the total amount of timesteps
    size = int((tEnd - t0) / dt + 1)

    # Calculate the amount of timesteps each processor is responsible for
    [chunks, offsets] = chunks_and_offsets(nProcs, size)

    # Perform the generation
    for i in range(chunks[rank]):
        t = t0 + dt * i + dt * int(offsets[rank])
        t = float(("{0:." + str(timePrecision) + "f}").format(t))
        position = int(offsets[rank]) + i

        if (rank == 0) and (np.mod(i, int(chunks[rank] / 10)) == 0):
            print("     Rescaled about " + str(int(i / chunks[rank] * 100)) + "%")

        # The random signal for all three components
        uTilde = np.random.normal(0, 1, (nPoints, 3))

        uX = uMean[:, 0] + np.sum(a[:, 0, :]*uTilde, axis=1)
        uY = uMean[:, 1] + np.sum(a[:, 1, :]*uTilde, axis=1)
        uW = uMean[:, 2] + np.sum(a[:, 2, :]*uTilde, axis=1)

        writerFunction(t, position, uX, uY, uW)

