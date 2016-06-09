import eddylicious
from eddylicious.writers.ofnative_writers import *
import numpy as np
from os import path


def test_point_writer(tmpdir):
    dsvDir = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                     "channel_flow_180", "dsv_output")

    pointsY = np.load(path.join(dsvDir, "pointsY.npy"))
    pointsZ = np.load(path.join(dsvDir, "pointsZ.npy"))
    xVal = 0.0
    writeFile = tmpdir.join("points")
    writePath = writeFile.strpath
    write_points_to_ofnative(writePath, pointsY, pointsZ, xVal)

    with file(writePath) as pointsFile:
        writtenPoints = [line.rstrip(')\n') for line in pointsFile]

    writtenPoints = [line.lstrip('(') for line in writtenPoints]
    writtenPoints = writtenPoints[9:-1]
    writtenPoints = np.genfromtxt(writtenPoints)

    assert np.all(writtenPoints[:, 0] ==
                  xVal*np.ones(writtenPoints[:, 0].shape))
    assert np.all(writtenPoints[:, 1] ==
                  pointsY.reshape((pointsY.size, -1), order='F')[:, 0])

    assert np.all(writtenPoints[:, 2] ==
                  pointsZ.reshape((pointsY.size, -1), order='F')[:, 0])


def test_velocity_writer(tmpdir):
    dsvDir = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180", "dsv_output")

    uX = np.load(path.join(dsvDir, "1000.01", "uX.npy"))
    uY = np.load(path.join(dsvDir, "1000.01", "uY.npy"))
    uZ = np.load(path.join(dsvDir, "1000.01", "uZ.npy"))

    writeFile = tmpdir.mkdir("u")
    writePath = writeFile.strpath
    write_velocity_to_ofnative(writePath, "0.1", uX, uY, uZ)

    with file(path.join(writePath, "0.1", "U")) as uFile:
        writtenU= [line.rstrip(')\n') for line in uFile]

    writtenU = [line.lstrip('(') for line in writtenU]
    writtenU = writtenU[10:-1]
    writtenU = np.genfromtxt(writtenU)

    uX = uX.reshape((uX.size, -1), order='F')[:, 0]
    uY = uY.reshape((uY.size, -1), order='F')[:, 0]
    uZ = uZ.reshape((uZ.size, -1), order='F')[:, 0]

    assert np.all(writtenU[:, 0] == uX)
    assert np.all(writtenU[:, 1] == uY)
    assert np.all(writtenU[:, 2] == uZ)
