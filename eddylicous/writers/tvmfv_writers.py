import os
import numpy as np

"""Functions for writing to the native format of the
timeVaryingMappedFixedValue boundary in OpenFOAM.
"""


def write_points_to_tvmfv(writePath, pointsY, pointsZ, xVal):
    """Write the points in a format used by OpenFOAM's
    timeVaryingMappedFixedValue boundary conidition.
    """

    pointsHeader = \
        "FoamFile\n{\nversion 2.0;\nformat ascii;\n\
        class vectorField;\nobject values;\n}\n"

    points = np.zeros((pointsY.size, 2))
    points[:, 0] = np.reshape(pointsY, (pointsY.size, -1), order='F')[:, 0]
    points[:, 1] = np.reshape(pointsZ, (pointsZ.size, -1), order='F')[:, 0]
    points = np.concatenate((xVal*np.ones((points.shape[0], 1)), points),
                            axis=1)
    np.savetxt(writePath, points,
               header=pointsHeader+str(points.shape[0])+"\n(", footer=")\n",
               comments="", fmt='(%f %f %f)')


def write_u_to_tvmfv(writePath, t, u):
    """Write the velocity field in a format used by OpenFOAM's
    timeVaryingMappedFixedValue boundary conidition.
    """

    vectorHeader = \
        "FoamFile\n{\nversion 2.0;\nformat ascii;\n\
        class vectorAverageField;\nobject points;\n}\n\n(0 0 0)\n"

    if not os.path.exists(os.path.join(writePath, str(t))):
        os.mkdir(os.path.join(writePath, str(t)))

    np.savetxt(os.path.join(writePath, str(t), "U"), u,
               header=vectorHeader+str(u.shape[0])+"\n(", footer=")\n",
               comments="", fmt='(%f %f %f)')
