"""Functions for writing to the native format of the
timeVaryingMappedFixedValue boundary in OpenFOAM.

"""

import numpy as np
import h5py as h5py


__all__ = ["write_points_to_hdf5", "write_u_to_hdf5"]



def write_points_to_hdf5(writePath, pointsY, pointsZ, xVal):
    """Write the points into a HDF5 file.

    This function will save the points into a HDF5 file.
    The points will be transformed into 1d arrays.
    The resulting dataset is called points and lies in the
    root of the file.

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


def write_u_to_hdf5(file, t, u, iter, size):
    """Write the velocity field into an HDF5 file.

    Will add datasets "time" and "velocity"

    Parameters
    ---------
    writePath : str
        The path of the HDF5 file.
    t : float
        The value of time associated with the written
        velocity field.
    u : ndarray
        Array containing the velocity field.
    iter: int
        The position of along the time axis.
    size: int
        The total size of the time axis.
    """

    if (iter > size-1):
        print "WARNING in write_u_to_hdf5. Write position larger \
              than total database size. Not writing."
        return

    file["time"][iter] = t
    file["velocity"][iter, :, :] = u

