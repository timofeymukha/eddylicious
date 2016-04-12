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
    [pY, pZ, yI, zI] = read_points_from_foamfile(readPath)

    assert np.all(load_points[1] == pY)
    assert np.all(load_points[2] == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_exclude_top(load_points):
    n = 10

    nPointsY = load_points[1].shape[0]
    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_points_from_foamfile(readPath,
                                                 excludeTop=nPointsY-n)

    assert np.all(load_points[1][:n, :] == pY)
    assert np.all(load_points[2][:n, :] == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)

def test_read_points_exclude_bot(load_points):
    n = 10

    nPointsY = load_points[1].shape[0]
    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    [pY, pZ, yI, zI] = read_points_from_foamfile(readPath,
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

    [pY, pZ, yI, zI] = read_points_from_foamfile(readPath, addValBot=0,
                                                 addValTop=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


def test_read_points_add_zeros_bot_exclude_top_midvalue(load_points):
    n = 10

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)[:n, :]
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2], axis=0)[:n, :]
    pointsY[-1, :] = 1.0

    readPath = path.join(load_points[0], "foam_file_output", "1000.01",
                         "faceCentres")

    nPointsY = load_points[1].shape[0]+1

    [pY, pZ, yI, zI] = read_points_from_foamfile(readPath, addValBot=0,
                                                 excludeTop=nPointsY-n,
                                                 midValue=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(load_points[3] == yI)
    assert np.all(load_points[4] == zI)


# Tests for the velocity reader
def test_read_velocity_add_zeros_bot_some_points_interpolate():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    uX = np.load(path.join(prefix, "dsv_output", "1000.03", "uX.npy"))
    uY = np.load(path.join(prefix, "dsv_output", "1000.03", "uY.npy"))
    uZ = np.load(path.join(prefix, "dsv_output", "1000.03", "uZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    nPointsY = 15

    uX = np.append(np.zeros((1, 72)), uX, axis=0)
    uY = np.append(np.zeros((1, 72)), uY, axis=0)
    uZ = np.append(np.zeros((1, 72)), uZ, axis=0)

    uX[nPointsY-1, :] = 0.5*(uX[nPointsY-2, :] + uX[nPointsY, :])
    uY[nPointsY-1, :] = 0.5*(uY[nPointsY-2, :] + uY[nPointsY, :])
    uZ[nPointsY-1, :] = 0.5*(uZ[nPointsY-2, :] + uZ[nPointsY, :])

    readFunc = read_velocity_from_foamfile(path.join(prefix,
                                                     "foam_file_output"), "",
                                           nPointsY,
                                           72, yInd, zInd, addValBot=0,
                                           interpolate=1)
    [uXR, uYR, uZR] = readFunc(1000.03)

    assert np.all(uX[:nPointsY, :] == uXR[:nPointsY, :])
    assert np.all(uY[:nPointsY, :] == uYR[:nPointsY, :])
    assert np.all(uZ[:nPointsY, :] == uZR[:nPointsY, :])


def test_read_velocity_add_zeros_bot_some_points_no_interpolate():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    uX = np.load(path.join(prefix, "dsv_output", "1000.01", "uX.npy"))
    uY = np.load(path.join(prefix, "dsv_output", "1000.01", "uY.npy"))
    uZ = np.load(path.join(prefix, "dsv_output", "1000.01", "uZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    nPointsY = 15

    uX = np.append(np.zeros((1, 72)), uX, axis=0)
    uY = np.append(np.zeros((1, 72)), uY, axis=0)
    uZ = np.append(np.zeros((1, 72)), uZ, axis=0)

    readFunc = read_velocity_from_foamfile(path.join(prefix,
                                                     "foam_file_output"), "",
                                           nPointsY, 72, yInd, zInd, addValBot=0,
                                           interpolate=0)

    [uXR, uYR, uZR] = readFunc(1000.01)

    assert np.all(uX[:nPointsY, :] == uXR[:nPointsY, :])
    assert np.all(uY[:nPointsY, :] == uYR[:nPointsY, :])
    assert np.all(uZ[:nPointsY, :] == uZR[:nPointsY, :])


def test_read_velocity_add_zeros_bot_top_all_points_no_interpolate():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    uX = np.load(path.join(prefix, "dsv_output", "1000.01", "uX.npy"))
    uY = np.load(path.join(prefix, "dsv_output", "1000.01", "uY.npy"))
    uZ = np.load(path.join(prefix, "dsv_output", "1000.01", "uZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))


    uX = np.append(np.zeros((1, 72)), uX, axis=0)
    uY = np.append(np.zeros((1, 72)), uY, axis=0)
    uZ = np.append(np.zeros((1, 72)), uZ, axis=0)

    uX = np.append(uX, np.zeros((1, 72)), axis=0)
    uY = np.append(uY, np.zeros((1, 72)), axis=0)
    uZ = np.append(uZ, np.zeros((1, 72)), axis=0)

    nPointsY = uX.shape[0]

    readFunc = read_velocity_from_foamfile(path.join(prefix,
                                                     "foam_file_output"), "",
                                           nPointsY, 72, yInd, zInd, addValBot=0,
                                           addValTop=0, interpolate=0)

    [uXR, uYR, uZR] = readFunc(1000.01)



    assert np.all(uX[:, :] == uXR[:, :])
    assert np.all(uY[:, :] == uYR[:, :])
    assert np.all(uZ[:, :] == uZR[:, :])




