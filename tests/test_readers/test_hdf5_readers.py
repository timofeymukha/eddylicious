import eddylicious
from eddylicious.readers.hdf5_readers import *
import numpy as np
from os import path
import os
import pytest
import h5py


@pytest.fixture
def load_points():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    pointsY = np.load(path.join(prefix, "dsv_output", "pointsY.npy"))
    pointsZ = np.load(path.join(prefix, "dsv_output", "pointsZ.npy"))
    return [prefix, pointsY, pointsZ]

@pytest.fixture
def create_hdf5(tmpdir, load_points):
    writeFile = tmpdir.join("test.hdf5")
    writePath = writeFile.strpath
    if path.isfile(writePath):
        os.remove(writePath)
    dbFile = h5py.File(writePath, 'a')
    group = dbFile.create_group("points")
    group.create_dataset("pointsY", data=load_points[1])
    group.create_dataset("pointsZ", data=load_points[2])
    return writePath


# Do not add zeros on top and bottom, take all points along y
def test_read_points_no_add_zeros_all_points(load_points, create_hdf5):


    [pY, pZ] = read_points_from_hdf5(create_hdf5)

    assert np.all(load_points[1] == pY)
    assert np.all(load_points[2] == pZ)


def test_read_points_exclude_top(load_points, create_hdf5):
    n = 10

    nPointsY = load_points[1].shape[0]

    [pY, pZ,] = read_points_from_hdf5(create_hdf5,
                                                 excludeTop=nPointsY-n)

    assert np.all(load_points[1][:n, :] == pY)
    assert np.all(load_points[2][:n, :] == pZ)

def test_read_points_exclude_bot(load_points, create_hdf5):
    n = 10

    [pY, pZ,] = read_points_from_hdf5(create_hdf5,
                                      excludeBot=n)

    assert np.all(load_points[1][n:, :] == pY)
    assert np.all(load_points[2][n:, :] == pZ)


def test_read_points_add_zeros_bot_top(load_points, create_hdf5):

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2],
                        axis=0)
    pointsY = np.append(pointsY, np.zeros((1, 72)), axis=0)
    pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

    [pY, pZ] = read_points_from_hdf5(create_hdf5, addValBot=0, addValTop=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)


def test_read_points_add_zeros_bot_exclude_top_midvalue(load_points,
                                                        create_hdf5):
    n = 10

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)[:n, :]
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2],
                        axis=0)[:n, :]
    pointsY[-1, :] = 1.0

    nPointsY = load_points[1].shape[0]+1

    [pY, pZ] = read_points_from_hdf5(create_hdf5, addValBot=0,
                                                 excludeTop=nPointsY-n,
                                                 midValue=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)

