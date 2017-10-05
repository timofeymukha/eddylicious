# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for reading fields stored in the hdf5 format.

"""
import numpy as np
import h5py

__all__ = ["read_structured_points_hdf5", "read_structured_velocity_hdf5"]


def read_structured_points_hdf5(readPath, addValBot=float('nan'),
                                addValTop=float('nan'), excludeBot=0,
                                excludeTop=0, exchangeValBot=float('nan'),
                                exchangeValTop=float('nan')):
    """Read the coordinates of the points from a hdf5 file.


    Reads in the locations of the face centers, stored in  a hdf5 file.

    The function supports manipulating the points in certain ways, see
    the parameter list below.


    Parameters
    ----------
    readPath : str
        The path to the file containing the points
    addValBot : float, optional
        Whether to append a row of zeros from below.
    addValTop : float, optional
        Whether to append a row of zeros from above.
    excludeBot : int, optional
        How many points to remove from the bottom in the y direction.
        (default 0).
    excludeTop: int, optional
        How many points to remove from the top in the y direction.
        (default 0).
    exchangeValBot : float, optional
        Exchange the value of y at the bottom.
    exchangeValTop : float, optional
        Exchange the value of y at the top.


    Returns
    -------
    List of ndarrays
        The list contains 2 items

        pointsY :
        A 2d ndarray containing the y coordinates of the points.

        pointsZ :
        A 2d ndarray containing the z coordinates of the points.

    """
    dbFile = h5py.File(readPath, 'a')

    pointsY = dbFile["points"]["pointsY"]
    pointsZ = dbFile["points"]["pointsZ"]

    nPointsZ = pointsY.shape[1]

    # Add points at y = 0 and y = max(y)
    if not np.isnan(addValBot):
        pointsY = np.append(addValBot*np.ones((1, nPointsZ)), pointsY, axis=0)
        pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)
    if not np.isnan(addValTop):
        pointsY = np.append(pointsY, addValTop*np.ones((1, nPointsZ)), axis=0)
        pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

    # Cap the points

    nPointsY = pointsY.shape[0]
    if excludeTop:
        pointsY = pointsY[:(nPointsY-excludeTop), :]
        pointsZ = pointsZ[:(nPointsY-excludeTop), :]

    if excludeBot:
        pointsY = pointsY[excludeBot:, :]
        pointsZ = pointsZ[excludeBot:, :]

    if not np.isnan(exchangeValBot):
        pointsY[0, :] = exchangeValBot

    if not np.isnan(exchangeValTop):
        pointsY[-1, :] = exchangeValTop

    return [pointsY, pointsZ]


def read_structured_velocity_hdf5(readPath,
                                  addValBot=(float('nan'), float('nan'),
                                                       float('nan')),
                                  addValTop=(float('nan'), float('nan'),
                                       float('nan')),
                                  excludeBot=0, excludeTop=0,
                                  interpValBot=False, interpValTop=False):
    """ Read the values of the velocity field from a foamFile-format
    file.

    Reads in the values of the velocity components stored as in hdf5
    file format.

    Some manipulation with the read-in data is also available via the
    optional parameters.

    Parameters
    ---------
    readPath : str
        The path to the file containing the velocity field.
    addValBot : tuple of three floats, optional
        Append a row of values from below.
    addValTop : tuple of three floats, optional
        Append a row of values from above.
    excludeBot : int, optional
        How many points to remove from the bottom in the y direction.
        (default 0).
    excludeTop: int, optional
        How many points to remove from the top in the y direction.
        (default 0).
    interpValBot : bool, optional
        Whether to interpolate the first value in the wall-normal
        direction using two points. (default False)
    interpValTop : bool, optional
        Whether to interpolate the last value in the wall-normal
        direction using two points. (default False)

    Returns
    -------
    function
        A function of one variable (the time-index) that will actually
        perform the reading.
    """

    def read(timeIndex):
        """
        A function that will actually perform the reading.

        Parameters
        ----------
        timeIndex: int
            The value of the time-index, i.e. the location in the
            times-array.

        Returns
        -------
        List of ndarrays
            The list contains three items, each a 2d array,
            corresponding to the three components of the velocity field,
            the order of the components in the list is x, y and the z.

        """
        dbFile = h5py.File(readPath, 'r')

        uX = dbFile["velocity"]["uX"][timeIndex, :, :]
        uY = dbFile["velocity"]["uY"][timeIndex, :, :]
        uZ = dbFile["velocity"]["uZ"][timeIndex, :, :]

        nPointsZ = uX.shape[1]

        # SAME VALUE FOR ALL COMPONENTS?! Change to list..
        if not np.isnan(addValBot[0]):
            uX = np.append(addValBot[0]*np.ones((1, nPointsZ)), uX, axis=0)
            uY = np.append(addValBot[1]*np.ones((1, nPointsZ)), uY, axis=0)
            uZ = np.append(addValBot[2]*np.ones((1, nPointsZ)), uZ, axis=0)

        if not np.isnan(addValTop[0]):
            uX = np.append(uX, addValTop[0]*np.ones((1, nPointsZ)), axis=0)
            uY = np.append(uY, addValTop[1]*np.ones((1, nPointsZ)), axis=0)
            uZ = np.append(uZ, addValTop[2]*np.ones((1, nPointsZ)), axis=0)

        nPointsY = uX.shape[0]
        topmostPoint = nPointsY-excludeTop

        # Interpolate for the last point in the wall-normal direction
        if interpValTop and excludeTop:
            uX[topmostPoint-1, :] = 0.5*(uX[topmostPoint-1, :] +
                                         uX[topmostPoint, :])
            uY[topmostPoint-1, :] = 0.5*(uY[topmostPoint-1, :] +
                                         uY[topmostPoint, :])
            uZ[topmostPoint-1, :] = 0.5*(uZ[topmostPoint-1, :] +
                                         uZ[topmostPoint, :])

        # Interpolate for the first point in the wall-normal direction
        if interpValBot and excludeBot:
            uX[excludeBot, :] = 0.5*(uX[excludeBot-1, :] +
                                     uX[excludeBot, :])
            uY[excludeBot, :] = 0.5*(uY[excludeBot-1, :] +
                                     uY[excludeBot, :])
            uZ[excludeBot, :] = 0.5*(uZ[excludeBot-1, :] +
                                     uZ[excludeBot, :])

        # Cap the points
        if excludeTop:
            uX = uX[:topmostPoint, :]
            uY = uY[:topmostPoint, :]
            uZ = uZ[:topmostPoint, :]

        if excludeBot:
            uX = uX[excludeBot:, :]
            uY = uY[excludeBot:, :]
            uZ = uZ[excludeBot:, :]

        return [uX, uY, uZ]

    read.reader = "hdf5"
    return read
