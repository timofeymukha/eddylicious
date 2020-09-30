# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for writing to the native format of the
timeVaryingMappedFixedValue boundary in OpenFOAM.

"""
import os
import numpy as np

__all__ = ["write_points_to_ofnative", "write_velocity_to_ofnative"]


def write_points_to_ofnative(writePath, pointsY, pointsZ, xVal):
    """Write the points in a format used by OpenFOAM's
    timeVaryingMappedFixedValue boundary condition.

    Parameters
    ----------
    writePath : str
        The path where to write the points file. Should commonly be
        constant/boundaryData/nameOfInletPatch.
    pointsY : ndarray
        A 2d array containing the values of y for the face centres.
    pointsZ : ndarray
        A 2d array containing the values of z for the face centres.
    xVal : float
        The x-location of the inflow plane.

    """

    points = np.zeros((pointsY.size, 2))
    points[:, 0] = np.reshape(pointsY, (pointsY.size, -1), order='F')[:, 0]
    points[:, 1] = np.reshape(pointsZ, (pointsZ.size, -1), order='F')[:, 0]
    points = np.concatenate((xVal*np.ones((points.shape[0], 1)), points),
                            axis=1)
    np.savetxt(writePath, points,
               header=str(points.shape[0])+"\n(", footer=")\n",
               comments="", fmt='(%e %e %e)')


def write_velocity_to_ofnative(writePath, t, position, uX, uY, uZ):
    """Write the velocity field in a format used by OpenFOAM's
    timeVaryingMappedFixedValue boundary condition.

    Parameters
    ----------
    writePath : str
        The path where to write the time directories containing the U
        files. Commonly constant/boundaryData/nameOfInletPatch.
    t : float
        The value of time associated with the written velocity field.
    position : int
        Not used here, present to comply with common writer interface.
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

    if not os.path.exists(os.path.join(writePath, str(t))):
        os.mkdir(os.path.join(writePath, str(t)))

    uX = np.reshape(uX, (uX.size, -1), order='F')
    uY = np.reshape(uY, (uY.size, -1), order='F')
    uZ = np.reshape(uZ, (uZ.size, -1), order='F')

    u = np.concatenate((uX, uY, uZ), axis=1)

    np.savetxt(os.path.join(writePath, str(t), "U"), u,
               header=str(u.shape[0])+"\n(", footer=")\n",
               comments="", fmt='(%e %e %e)')
