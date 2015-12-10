import numpy as np

"""Functions for reading fields stored in the foamFile format"""


def read_points_from_foamfile(readPath, addZeros=1, nPointsY=0, delta=1):
    """Read the coordinates of the points from foamFile-format file.


    Reads in the locations of the face centers, stored in foamFile
    format by OpenFOAM, and transforms them into 2d numpy arrays.

    The locations are sorted along y and z.

    There is a possibility to add a coordinates for y=0.
    Also, one can keep only a part of the points in the y
    direction and exchange the last value of y to some
    delta.
    This is used when reading in precursor data from channel
    flow, since we are only interested in half of the channel.

    Returns a list with 4 items: pointsY, pointsZ, indY and indZ.
    The first two are 2d numpy array with the locations.
    The last two are the sorting indices.

    Keyword argments:
    addZeros -- whether to add coordinates for y=0 (default 1)
    nPointsY -- how many points to keep in the y direction.
               Zeros means all points are kept (default 0).
    delta -- the  value of the channel-half width. (default 1)
    """

    with file(readPath) as pointsFile:
        points = [line.rstrip(')\n') for line in pointsFile]

    points = [line.lstrip('(') for line in points]
    points = points[3:-1]
    points = np.genfromtxt(points)[:, 1:]

# Sort the points
# Sort along y first
    yInd = np.argsort(points[:, 0])
    points[:, 0] = points[yInd, 0]
    points[:, 1] = points[yInd, 1]

# Find the number of points along z
    nPointsZ = 0
    for i in xrange(points[:, 0].size):
        if points[i, 0] == points[0, 0]:
            nPointsZ += 1
        else:
            break

# Reshape into a 2d array
    pointsY = np.copy(np.reshape(points[:, 0], (-1, nPointsZ)))
    pointsZ = np.copy(np.reshape(points[:, 1], (-1, nPointsZ)))

# For each y order the points in z

    zInd = np.zeros(pointsZ.shape, dtype=np.int)

    for i in xrange(pointsZ.shape[0]):
        zInd[i, :] = np.argsort(pointsZ[i, :])
        pointsZ[i, :] = pointsZ[i, zInd[i, :]]


# Add points at y = 0
    if addZeros:
        pointsY = np.append(np.zeros((1, pointsY.shape[1])), pointsY, axis=0)
        pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)

# Cap the points, to include ony one half of the channel
# Makes y=delta the last point
    if nPointsY:
        pointsY = pointsY[:nPointsY, :]
        pointsZ = pointsZ[:nPointsY, :]
        pointsY[-1, :] = delta
    else:
        nPointsY = pointsY.shape[0]

    return [pointsY, pointsZ, yInd, zInd]


def read_u_from_foamfile(readPath, nPointsY, nPointsZ, yInd, zInd):
    """ Read the values of the velocity field from a foamFile-format file.

        Reads in the values of the velocity components stored
        as in foamFile file-format.
        Requires the number of points to keep in the y
        direction and the sorting indices, so that one can
        sort in the same way as the face-centers are sorted.

        Returns a list with three items: numpy arrays for the
        three components of velocity.
    """

    with file(readPath) as UFile:
        u = [line.rstrip(')\n') for line in UFile]

    u = [line.lstrip('(') for line in u]
    u = u[3:-1]
    u = np.genfromtxt(u)

    # Sort along y
    u[:, 0] = u[yInd, 0]
    u[:, 1] = u[yInd, 1]
    u[:, 2] = u[yInd, 2]

    # Reshape to 2d
    uX = np.copy(np.reshape(u[:, 0], (-1, nPointsZ)))
    uY = np.copy(np.reshape(u[:, 1], (-1, nPointsZ)))
    uZ = np.copy(np.reshape(u[:, 2], (-1, nPointsZ)))

    # Sort along z
    for i in xrange(uX.shape[0]):
        uX[i, :] = uX[i, zInd[i, :]]
        uY[i, :] = uY[i, zInd[i, :]]
        uZ[i, :] = uZ[i, zInd[i, :]]

    # Prepend with zeros
    uX = np.append(np.zeros((1, nPointsZ)), uX, axis=0)
    uY = np.append(np.zeros((1, nPointsZ)), uY, axis=0)
    uZ = np.append(np.zeros((1, nPointsZ)), uZ, axis=0)

    # Interpolate to get data at y=delta
    uX[nPointsY-1, :] = 0.5*(uX[nPointsY-2, :] + uX[nPointsY, :])
    uY[nPointsY-1, :] = 0.5*(uY[nPointsY-2, :] + uY[nPointsY, :])
    uZ[nPointsY-1, :] = 0.5*(uZ[nPointsY-2, :] + uZ[nPointsY, :])

    # Remove data above y=delta
    uX = uX[:nPointsY, :]
    uY = uY[:nPointsY, :]
    uZ = uZ[:nPointsY, :]

    return [uX, uY, uZ]
