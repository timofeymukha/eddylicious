import numpy as np
import h5py as h5py

"""Functions for reading fields stored in the hdf5 format"""


def read_points_from_hdf5(readPath, addZeros=0, nPointsY=0, delta=1):
    """Read the coordinates of the points from a hdf5 file.


    Reads in the locations of the face centers, stored
    in  a hdf5 file.

    The function supports considering only a number of
    points in the wall-normal direction and exchanging
    the last wall-normal position with the value of the
    half-width of the channel. 
    Also, adding a row of zeros as the first wall-normal
    position is possible.

    This is convenient when rescaling from channel flow
    is performed using Lund et al's method, which requires
    a liner interpolant across the domain half-width.
    Adding the value at the center of the channel and at
    the wall, which are otherwise absent on a finite volume
    grid, insures that the interpolant will cover the whole
    interval [0, delta].


    Parameters
    ----------
    readPath : str
        The path to the file containing the points
    addZeros : bool, optional,
        Whether to add coordinates for y=0 (default is 1)
    nPointsY : int, optional
        How many points to keep in the y direction. Zero
        means all points are kept (default 0).
        
        If nPointsY is provided, the last value in the wall-
        normal direction is exchanged to value of delta,
        see below.
    delta : bool, optional
        The value of the channel-half width. (default 1).
        Must be provided if nPoints is provided.


    Returns
    -------
    List of ndarrays
        The list contains 4 items
        pointsY :
            A 2d ndarray containing the y coordinates of the
            points.
        pointsZ :
            A 2d ndarray containing the z coordinates of the
            points.
    """

    dbFile = h5py.File(readPath, 'a')

    pointsY = dbFile["points"]["pointsY"]
    pointsZ = dbFile["points"]["pointsZ"]

    nPointsZ  = pointsY.shape[1]

# Add points at y = 0 and y = max(y)
    if addZeros:
        pointsY = np.append(np.zeros((1, nPointsZ)), pointsY, axis=0)
        pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)
        if nPointsY == 0:
            pointsY = np.append(pointsY, np.zeros((1, nPointsZ)), axis=0)
            pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

# Cap the points, to include ony one half of the channel
# Makes y=delta the last point
    if nPointsY:
        pointsY = pointsY[:nPointsY, :]
        pointsZ = pointsZ[:nPointsY, :]
        pointsY[-1, :] = delta
    else:
        nPointsY = pointsY.shape[0]

    return [pointsY, pointsZ]


def read_u_from_hdf5(readPath, timeIndex):
    """ Read the values of the velocity field from a foamFile-format file.

    Reads in the values of the velocity components stored
    as in hdf5 file format.

    Parameters
    ---------
    readPath : str
        The path to the file containing the velocity field.
    timeIndex : int
        The time-index associated with the data, i.e. the index of the
        first dimension of the hdf5 file.


    Returns
    -------
    List of ndarrays
        The list contains three items, each a 2d array,
        corresponding to the three components of the
        velocity field, the order of the components in the
        list is x, y and the z.
    """
    dbFile = h5py.File(readPath, 'a')

    uX = dbFile["velocity"]["uX"][timeIndex, :, :]
    uY = dbFile["velocity"]["uY"][timeIndex, :, :]
    uZ = dbFile["velocity"]["uZ"][timeIndex, :, :]

    # Prepend with zeros
    #uX = np.append(np.zeros((1, nPointsZ)), uX, axis=0)
    #uY = np.append(np.zeros((1, nPointsZ)), uY, axis=0)
    #uZ = np.append(np.zeros((1, nPointsZ)), uZ, axis=0)

    # Interpolate to get data at y=delta
    #uX[nPointsY-1, :] = 0.5*(uX[nPointsY-2, :] + uX[nPointsY, :])
    #uY[nPointsY-1, :] = 0.5*(uY[nPointsY-2, :] + uY[nPointsY, :])
    #uZ[nPointsY-1, :] = 0.5*(uZ[nPointsY-2, :] + uZ[nPointsY, :])

    # Remove data above y=delta
    #uX = uX[:nPointsY, :]
    #uY = uY[:nPointsY, :]
    #uZ = uZ[:nPointsY, :]

    return [uX, uY, uZ]
