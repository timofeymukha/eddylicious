import numpy as np
import h5py as h5py

"""Functions for writing to the native format of the
timeVaryingMappedFixedValue boundary in OpenFOAM.
"""


def write_points_to_hdf5(writePath, pointsY, pointsZ, xVal):
    """Write the points into a HDF5 file.

    This function will save the points into a HDF5 file.
    The points will be transformed into 1d arrays.
    The resulting dataset is called points and lies in the
    root of the file.

    Parameters
    ----------
    writePath : str The path of the HDF5 file.
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


def write_u_to_hdf5(writePath, t, u):
    """Write the velocity field into an HDF5 file.

    Parameters
    ---------
    writePath : str
        The path of the HDF5 file.
    t : float
        The value of time associated with the written
        velocity field.
    u : ndarray
        Array containing the velocity field
    """

    dbFile = h5py.File(writePath, 'a')

    time = "/time" in dbFile

    if time:
        size = dbFile["time"].size
        dbFile["time"].resize(size+1, axis=0)
        dbFile["time"][-1] = t
    else:
        time = dbFile.create_dataset("time", data=t*np.ones((1, 1)),
                                     chunks=True, maxshape=(None, 1))

    velocity = "/velocity" in dbFile

    if velocity:
        size = dbFile["velocity"].shape[2]
        dbFile["velocity"].resize(size+1, axis=2)
        dbFile["velocity"][:,:,-1] = u
    else:
        velocity = dbFile.create_dataset(
            "velocity",
            data=u[:, :, np.newaxis]*np.ones((u.shape[0], u.shape[1], 1)),
            chunks=True, maxshape=(u.shape[0], u.shape[1], None))

    dbFile.close()
