import eddylicious
from eddylicious.readers.foamfile_readers import *
import numpy as np
from os import path
import pytest


@pytest.fixture
def load_points(scope="module"):
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))
    return [prefix, pointsY, pointsZ, yInd, zInd]


def test_read_points_all_points(load_points):

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")
    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath)

    assert np.all(load_points[1] == pY)
    assert np.all(load_points[2] == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_exclude_top(load_points):
    n = 10

    nPointsY = load_points[1].shape[0]
    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath,
                                                       excludeTop=nPointsY-n)

    assert np.all(load_points[1][:n, :] == pY)
    assert np.all(load_points[2][:n, :] == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_exclude_bot(load_points):
    n = 10

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath,
                                                       excludeBot=n)

    assert np.all(load_points[1][n:, :] == pY)
    assert np.all(load_points[2][n:, :] == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_add_zeros_bot_top(load_points):

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2],
                        axis=0)
    pointsY = np.append(pointsY, np.zeros((1, 72)), axis=0)
    pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath, addValBot=0,
                                                       addValTop=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_add_zeros_bot_exclude_top_interp_top(load_points):
    n = 10

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)[:n, :]
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2], axis=0)[:n, :]
    pointsY[-1, :] = 1.0

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    nPointsY = load_points[1].shape[0]+1

    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath, addValBot=0,
                                                       excludeTop=nPointsY-n,
                                                       exchangeValTop=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_add_zeros_top_exclude_bot_interp_bot(load_points):
    n = 10

    pointsY = np.append(load_points[1], np.zeros((1, 72)), axis=0)[n:, :]
    pointsZ = np.append(load_points[2],
                        np.array([load_points[2][0, :]]),
                        axis=0)[n:, :]
    pointsY[0, :] = 1.0

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_structured_points_foamfile(readPath, addValTop=0,
                                                       excludeBot=n,
                                                       exchangeValBot=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


# Tests for the velocity reader
@pytest.fixture
def load_vel(scope="module"):
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    uX = np.load(path.join(prefix, "dsv_output", "1000.01", "uX.npy"))
    uY = np.load(path.join(prefix, "dsv_output", "1000.01", "uY.npy"))
    uZ = np.load(path.join(prefix, "dsv_output", "1000.01", "uZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))
    return [prefix, uX, uY, uZ, yInd, zInd]


def test_read_velocity(load_vel):

    readFunc = read_structured_velocity_foamfile(path.join(load_vel[0],
                                                     "foam_file_output"), "",
                                                 72, load_vel[4], load_vel[5])
    [uXR, uYR, uZR] = readFunc(1000.01)

    assert np.all(load_vel[1] == uXR)
    assert np.all(load_vel[2] == uYR)
    assert np.all(load_vel[3] == uZR)


def test_read_velocity_add_zeros_bot_exclude_top_interp_top(load_vel):

    nPointsY = 15

    uX = np.append(np.zeros((1, 72)), load_vel[1], axis=0)
    uY = np.append(np.zeros((1, 72)), load_vel[2], axis=0)
    uZ = np.append(np.zeros((1, 72)), load_vel[3], axis=0)

    uX[nPointsY-1, :] = 0.5*(uX[nPointsY-1, :] + uX[nPointsY, :])
    uY[nPointsY-1, :] = 0.5*(uY[nPointsY-1, :] + uY[nPointsY, :])
    uZ[nPointsY-1, :] = 0.5*(uZ[nPointsY-1, :] + uZ[nPointsY, :])

    uX = uX[:nPointsY, :]
    uY = uY[:nPointsY, :]
    uZ = uZ[:nPointsY, :]

    readFunc = \
        read_structured_velocity_foamfile(path.join(load_vel[0],
                                              "foam_file_output"), "",
                                          72, load_vel[4], load_vel[5],
                                          addValBot=(0, 0, 0),
                                          excludeTop=load_vel[1].shape[0]-nPointsY+1,
                                          interpValTop=1)
    [uXR, uYR, uZR] = readFunc(1000.01)

    assert np.all(uX == uXR)
    assert np.all(uY == uYR)
    assert np.all(uZ == uZR)


def test_read_velocity_add_zeros_top_exclude_bot_interp_bot(load_vel):

    nPointsY = 15

    uX = np.append(load_vel[1], np.zeros((1, 72)), axis=0)
    uY = np.append(load_vel[2], np.zeros((1, 72)), axis=0)
    uZ = np.append(load_vel[3], np.zeros((1, 72)), axis=0)

    uX[nPointsY, :] = 0.5*(uX[nPointsY-1, :] + uX[nPointsY, :])
    uY[nPointsY, :] = 0.5*(uY[nPointsY-1, :] + uY[nPointsY, :])
    uZ[nPointsY, :] = 0.5*(uZ[nPointsY-1, :] + uZ[nPointsY, :])

    uX = uX[nPointsY:, :]
    uY = uY[nPointsY:, :]
    uZ = uZ[nPointsY:, :]

    readFunc = read_structured_velocity_foamfile(path.join(load_vel[0],
                                                     "foam_file_output"), "",
                                                 72, load_vel[4], load_vel[5],
                                                 addValTop=(0, 0 , 0),
                                                 excludeBot=nPointsY,
                                                 interpValBot=1)
    [uXR, uYR, uZR] = readFunc(1000.01)

    assert np.all(uX == uXR)
    assert np.all(uY == uYR)
    assert np.all(uZ == uZR)


