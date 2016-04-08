"""Functions for reading fields stored in the foamFile format

"""
import numpy as np

__all__ = ["read_points_from_foamfile", "read_u_from_foamfile"]


def read_points_from_foamfile(readPath, addZerosBot, addZerosTop,
                              nPointsY=0, midValue=0):
    """Read the coordinates of the points from a foamFile-format file.


    Reads in the locations of the face centers, stored in foamFile
    format by OpenFOAM, and transforms them into 2d numpy arrays.

    The points are sorted so that the axes of the arrays correspond to
    the wall-normal and spanwise directions. The points are first sorted
    along y, then reshaped into a 2d array and then sorted along z for
    each value of z.

    The function supports considering only a number of points in the
    wall-normal direction and exchanging the last wall-normal position
    with the value of the half-width of the channel. Also, adding a row
    of zeros as the first wall-normal position is possible.

    This is convenient when rescaling from channel flow is performed
    using Lund et al's method, which requires a liner interpolant across
    the domain half-width. Adding the value at the center of the channel
    and at the wall, which are otherwise absent on a finite volume grid,
    insures that the interpolant will cover the whole interval [0,
    delta].

    Parameters
    ----------
    readPath : str
        The path to the file containing the points.
    addZerosBot : bool
        Whether to append a row of zeros from below.
    addZerosTop : bool
        Whether to append a row of zeros from above.
    nPointsY : int, optional
        How many points to keep in the y direction. Zero means all
        points are kept (default 0).
    midValue : float, optional
        The value of the channel-half width. (default 0).

    Returns
    -------
    List of ndarrays
        The list contains 4 items
        pointsY :
        A 2d ndarray containing the y coordinates of the points.

        pointsZ :
        A 2d ndarray containing the z coordinates of the points.

        indY :
        The sorting indices from the sorting performed.

        indZ :
        The sorting indices from the sorting performed.

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


# Add points at y = 0 and y = max(y)
    if addZerosBot:
        pointsY = np.append(np.zeros((1, nPointsZ)), pointsY, axis=0)
        pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)
    if addZerosTop:
        pointsY = np.append(pointsY, np.zeros((1, nPointsZ)), axis=0)
        pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

# Cap the points, according to nPointsY
# Makes y=delta the last point
# NOTE! nPoints includes the added zeros
    if nPointsY:
        assert pointsY.shape[0] >= nPointsY

        pointsY = pointsY[:nPointsY, :]
        pointsZ = pointsZ[:nPointsY, :]

    if midValue:
        pointsY[-1, :] = midValue

    return [pointsY, pointsZ, yInd, zInd]


def read_u_from_foamfile(readPath, nPointsY, nPointsZ, yInd, zInd,
                         addZerosBot, addZerosTop, interpolate=0):
    """ Read the values of the velocity field from a foamFile-format
     file.

    Reads in the values of the velocity components stored as in foamFile
    file-format. The velocity field is read and the transformed into a
    2d numpy array, where the array's axes correspond to wall-normal and
    spanwise directions. To achieve this, the sorting indices obtained
    when reordering the mesh points are used.

    Parameters
    ----------
    interpolate
    readPath : str
        The path to the file containing the velocity field.
    nPointsY : int
        The number of points in the wall-normal direction to consider.
        For a single TBL it is enough to only consider half of the
        channel, for instance.
    nPointsZ : int
        The amount of points in the mesh in the spanwise direction.
    yInd : ndarray
        The sorting indices for sorting in the wall-normal direction.
    zInd : ndarray
        The sorting indices for sorting in the spanwise direction.
    addZerosBot : bool
        Whether to append a row of zeros from below.
    addZerosTop : bool
        Whether to append a row of zeros from above.
    interpolate : bool, optional
        Whether to interpolate the last value in the wall-normal
        direction using two points. Useful to get the center-value of
        the velocity from the channel flow.
    Returns
    -------
    List of ndarrays
        The list contains three items, each a 2d array,
        corresponding to the three components of the
        velocity field, the order of the components in the
        list is x, y and the z.

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

    if addZerosBot:
        uX = np.append(np.zeros((1, nPointsZ)), uX, axis=0)
        uY = np.append(np.zeros((1, nPointsZ)), uY, axis=0)
        uZ = np.append(np.zeros((1, nPointsZ)), uZ, axis=0)

    if addZerosTop:
        uX = np.append(uX, np.zeros((1, nPointsZ)), axis=0)
        uY = np.append(uY, np.zeros((1, nPointsZ)), axis=0)
        uZ = np.append(uZ, np.zeros((1, nPointsZ)), axis=0)

    # Interpolate for the last point in the wall-normal direction
    if interpolate:
        assert uX.shape[0] > nPointsY
        uX[nPointsY-1, :] = 0.5*(uX[nPointsY-2, :] + uX[nPointsY, :])
        uY[nPointsY-1, :] = 0.5*(uY[nPointsY-2, :] + uY[nPointsY, :])
        uZ[nPointsY-1, :] = 0.5*(uZ[nPointsY-2, :] + uZ[nPointsY, :])

    # Remove data above y=delta
    uX = uX[:nPointsY, :]
    uY = uY[:nPointsY, :]
    uZ = uZ[:nPointsY, :]

    return [uX, uY, uZ]
