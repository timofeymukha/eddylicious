# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for writing to an hdf5 file.

"""
import numpy as np
import h5py as h5py


__all__ = ["write_points_to_hdf5", "write_velocity_to_hdf5"]


def write_points_to_hdf5(hdf5File, pointsY, pointsZ, xVal):
    """Write the points into a HDF5 file.

    Savs the points into a HDF5 file. The points will be transformed
    into 1d arrays. The resulting dataset is called points and lies in
    the root of the file.

    Parameters
    ----------
    hdf5File : h5py.File
        The path of the HDF5 file.
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

    if "/points" in hdf5File:
        hdf5File.__delitem__("points")

    hdf5File.create_dataset("points", data=points)


def write_velocity_to_hdf5(hdf5File, t, uX, uY, uZ, position):
    """Write the velocity field into an HDF5 file.

    Will also write the corresponding time value.

    Parameters
    ----------
    hdf5File : h5py.File
        The the HDF5 file.
    t : float
        The value of time associated with the written
        velocity field.
    position: int
        The position of the data along the time axis.
    uX : ndarray
        A 2d ndarray containing the streamwise component of the velocity
        field.
    uY : ndarray
        A 2d ndarray containing the wall-normal component of the
        velocity
        field.
    uZ : ndarray
        A 2d ndarray containing the spanwise component of the velocity
        field.

    """

    uX = np.reshape(uX, (uX.size, -1), order='F')
    uY = np.reshape(uY, (uY.size, -1), order='F')
    uZ = np.reshape(uZ, (uZ.size, -1), order='F')

    u = np.concatenate((uX, uY, uZ), axis=1)

    size = hdf5File["time"].size

    if position >= size:
        raise ValueError("Write position larger than total database size.")

    hdf5File["time"][position] = t
    hdf5File["velocity"][position, :, :] = u
