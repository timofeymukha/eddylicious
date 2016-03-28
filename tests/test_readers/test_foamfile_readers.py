import eddylicious
from eddylicious.readers.foamfile_readers import *
import numpy as np
import pytest
from os import path
from numpy.testing import assert_almost_equal


# Do not add zeros on top and bottom, take all points along y
def test_read_points_no_add_zeros_all_points():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                     "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    [pY, pZ, yI, zI] = read_points_from_foamfile(path.join(prefix,
                                                           "foam_file_output",
                                                           "1000.01",
                                                           "inletSurface",
                                                           "faceCentres"),
                                                 addZeros=False, nPointsY=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(yInd == yI)
    assert np.all(zInd == zI)

# Do not add zeros on top and bottom, take part of points along y
def test_read_points_no_add_zeros_some_points():
    n = 10
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))[:n, :]
    pointsY[-1, :] = 1.0
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))[:n, :]
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    [pY, pZ, yI, zI] = read_points_from_foamfile(path.join(prefix,
                                                           "foam_file_output",
                                                           "1000.01",
                                                           "inletSurface",
                                                           "faceCentres"),
                                                 addZeros=False, nPointsY=n)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(yInd == yI)
    assert np.all(zInd == zI)

# Add zeros on top and bottom, take all points along y
def test_read_points_add_zeros_all_points():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    pointsY = np.append(np.zeros((1, 72)), pointsY, axis=0)
    pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)
    pointsY = np.append(pointsY, np.zeros((1, 72)), axis=0)
    pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

    [pY, pZ, yI, zI] = read_points_from_foamfile(path.join(prefix,
                                                           "foam_file_output",
                                                           "1000.01",
                                                           "inletSurface",
                                                           "faceCentres"),
                                                 addZeros=True, nPointsY=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(yInd == yI)
    assert np.all(zInd == zI)

# Add zeros at the bottom, take part of the along y
def test_read_points_add_zeros_some_points():
    n = 10
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))
    yInd = np.load(path.join(prefix, "dsv_output", "yInd.npy"))
    zInd = np.load(path.join(prefix, "dsv_output", "zInd.npy"))

    pointsY = np.append(np.zeros((1, 72)), pointsY, axis=0)[:n, :]
    pointsZ = np.append(np.array([pointsZ[0, :]]), pointsZ, axis=0)[:n, :]
    pointsY[-1, :] = 1.0


    [pY, pZ, yI, zI] = read_points_from_foamfile(path.join(prefix,
                                                           "foam_file_output",
                                                           "1000.01",
                                                           "inletSurface",
                                                           "faceCentres"),
                                                 addZeros=True, nPointsY=n)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)
    assert np.all(yInd == yI)
    assert np.all(zInd == zI)

