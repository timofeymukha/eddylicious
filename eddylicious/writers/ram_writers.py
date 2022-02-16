# This file is part of eddylicious
# (c) 2020 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for writing into the computer's memory, i.e. some variable

"""
import os
import numpy as np
from os.path import join

__all__ = ["write_data_to_ram"]


def write_data_to_ram(storage, t, position, uX, uY, uZ):
    """Write the velocity field to an ndarray stored in memory

    Useful for studying different generation methods rapidly, avoiding
    hard disk IO.

    Parameters
    ----------
    storage : ndarray
        The variable that will be stored 
    t : float
        The value of time associated with the written velocity field.
        Not used here, included to comply with the writer interface.
    position: int
        The position of the data along the time axis.
    uX : ndarray
        A 2d ndarray containing the streamwise component of the velocity
        field.
    uY : ndarray
        A 2d ndarray containing the wall-normal component of the
        velocity field.
    uZ : ndarray
        A 2d ndarray containing the spanwise component of the velocity
        field.

    """

    if not os.path.exists(writePath):
        os.mkdir(os.path.join(writePath))

    uX = np.reshape(uX, (uX.size, -1), order='F')
    uY = np.reshape(uY, (uY.size, -1), order='F')
    uZ = np.reshape(uZ, (uZ.size, -1), order='F')

    data = {"uX": uX[:, 0], "uY": uY[:, 0], "uZ": uZ[:, 0]}

    # Not sure why this is neeed, the points should already be ndarrays
    pointsY = np.array(pointsY)
    pointsZ = np.array(pointsZ)

    pointsX = xOrigin*np.ones(pointsY.size)

    fileName = join(writePath, "inflow_" + str(position))
    pyevtk.hl.pointsToVTK(fileName, pointsX, pointsY, pointsZ, data=data)
