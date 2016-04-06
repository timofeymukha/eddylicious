"""Functions for writing to the and hdf5 file.

"""
import numpy as np
import h5py as h5py


__all__ = ["write_points_to_hdf5", "write_velocity_to_hdf5"]


def write_points_to_hdf5(writePath, pointsY, pointsZ, xVal):
    """Write the points into a HDF5 file.

    This function will save the points into a HDF5 file. The points will
    be transformed into 1d arrays. The resulting dataset is called
    points and lies in the root of the file.

    Parameters
    ----------
    writePath : str
        The path of the HDF5 file.
    pointsY : ndarray
        A 2d array containing the values of y for the face
        centres.
    pointsZ : ndarray
        A 2d array containing the values of z for the face
        centres.
    xVal : float
        The x-location of the inflow plane.

    """
    dbFile = h5py.File(writePath, 'a')

    points = np.zeros((pointsY.size, 2))
    points[:, 0] = np.reshape(pointsY, (pointsY.size, -1), order='F')[:, 0]
    points[:, 1] = np.reshape(pointsZ, (pointsZ.size, -1), order='F')[:, 0]
    points = np.concatenate((xVal*np.ones((points.shape[0], 1)), points),
                            axis=1)

    if "/points" in dbFile:
        dbFile.__delitem__("points")

    dbFile.create_dataset("points", data=points)
    dbFile.close()


def write_velocity_to_hdf5(file, t, uX, uY, uZ, iteration):
    """Write the velocity field into an HDF5 file.

    Will add datasets "time" and "velocity"

    Parameters
    ---------
    writePath : str
        The path of the HDF5 file.
    t : float
        The value of time associated with the written
        velocity field.
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
    iteration: int
        The position of along the time axis.
    """

    uX = np.reshape(uX, (uX.size, -1), order='F')
    uY = np.reshape(uY, (uY.size, -1), order='F')
    uZ = np.reshape(uZ, (uZ.size, -1), order='F')

    u = np.concatenate((uX, uY, uZ), axis=1)

    size = file["time"].size

    if iteration >= size:
        raise ValueError("Write position larger than total database size.")

    file["time"][iteration] = t
    file["velocity"][iteration, :, :] = u

