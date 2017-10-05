# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

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
def load_vel():
    prefix = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180")
    uX = np.load(path.join(prefix, "dsv_output", "1000.01", "uX.npy"))
    uY = np.load(path.join(prefix, "dsv_output", "1000.01", "uY.npy"))
    uZ = np.load(path.join(prefix, "dsv_output", "1000.01", "uZ.npy"))
    return [prefix, uX, uY, uZ]

@pytest.fixture
def create_hdf5(tmpdir, load_points, load_vel):
    writeFile = tmpdir.join("test.hdf5")
    writePath = writeFile.strpath
    if path.isfile(writePath):
        os.remove(writePath)
    dbFile = h5py.File(writePath, 'a')
    group = dbFile.create_group("points")
    group.create_dataset("pointsY", data=load_points[1])
    group.create_dataset("pointsZ", data=load_points[2])
    groupVel = dbFile.create_group("velocity")
    groupVel.create_dataset("uX", (1, load_vel[1].shape[0],
                                   load_vel[1].shape[1]), dtype=np.float64)
    groupVel.create_dataset("uY", (1, load_vel[1].shape[0],
                                   load_vel[1].shape[1]), dtype=np.float64)
    groupVel.create_dataset("uZ", (1, load_vel[1].shape[0],
                                   load_vel[1].shape[1]), dtype=np.float64)

    dbFile["velocity"]["uX"][0, :, :] = load_vel[1]
    dbFile["velocity"]["uY"][0, :, :] = load_vel[2]
    dbFile["velocity"]["uZ"][0, :, :] = load_vel[3]

    return writePath


# Do not add zeros on top and bottom, take all points along y
def test_read_points_no_add_zeros_all_points(load_points, create_hdf5):


    [pY, pZ] = read_structured_points_hdf5(create_hdf5)

    assert np.all(load_points[1] == pY)
    assert np.all(load_points[2] == pZ)


def test_read_points_exclude_top(load_points, create_hdf5):
    n = 10

    nPointsY = load_points[1].shape[0]

    [pY, pZ,] = read_structured_points_hdf5(create_hdf5,
                                            excludeTop=nPointsY-n)

    assert np.all(load_points[1][:n, :] == pY)
    assert np.all(load_points[2][:n, :] == pZ)

def test_read_points_exclude_bot(load_points, create_hdf5):
    n = 10

    [pY, pZ,] = read_structured_points_hdf5(create_hdf5,
                                            excludeBot=n)

    assert np.all(load_points[1][n:, :] == pY)
    assert np.all(load_points[2][n:, :] == pZ)


def test_read_points_add_zeros_bot_top(load_points, create_hdf5):

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2],
                        axis=0)
    pointsY = np.append(pointsY, np.zeros((1, 72)), axis=0)
    pointsZ = np.append(pointsZ, np.array([pointsZ[0, :]]),  axis=0)

    [pY, pZ] = read_structured_points_hdf5(create_hdf5, addValBot=0, addValTop=0)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)


def test_read_points_add_zeros_bot_exclude_top_interp_top(load_points,
                                                          create_hdf5):
    n = 10

    pointsY = np.append(np.zeros((1, 72)), load_points[1], axis=0)[:n, :]
    pointsZ = np.append(np.array([load_points[2][0, :]]), load_points[2],
                        axis=0)[:n, :]
    pointsY[-1, :] = 1.0

    nPointsY = load_points[1].shape[0]+1

    [pY, pZ] = read_structured_points_hdf5(create_hdf5, addValBot=0,
                                           excludeTop=nPointsY-n,
                                           exchangeValTop=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)

def test_read_points_add_zeros_top_exclude_bot_interp_bot(load_points,
                                                          create_hdf5):
    n = 10

    pointsY = np.append(load_points[1], np.zeros((1, 72)), axis=0)[n:, :]
    pointsZ = np.append(load_points[2],
                        np.array([load_points[2][0, :]]),
                        axis=0)[n:, :]
    pointsY[0, :] = 1.0

    [pY, pZ] = read_structured_points_hdf5(create_hdf5, addValTop=0,
                                           excludeBot=n,
                                           exchangeValBot=1)

    assert np.all(pointsY == pY)
    assert np.all(pointsZ == pZ)


# Tests for the velocity reader
def test_read_velocity(load_vel, create_hdf5):

    readFunc = read_structured_velocity_hdf5(create_hdf5)

    [uXR, uYR, uZR] = readFunc(0)

    assert np.all(load_vel[1] == uXR)
    assert np.all(load_vel[2] == uYR)
    assert np.all(load_vel[3] == uZR)


def test_read_velocity_add_zeros_bot_exclude_top_interp_top(load_vel,
                                                            create_hdf5):
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
        read_structured_velocity_hdf5(create_hdf5, addValBot=(0, 0, 0),
                                      excludeTop=load_vel[1].shape[0]-nPointsY+1,
                                      interpValTop=1)
    [uXR, uYR, uZR] = readFunc(0)

    assert np.all(uX == uXR)
    assert np.all(uY == uYR)
    assert np.all(uZ == uZR)


def test_read_velocity_add_zeros_top_exclude_bot_interp_bot(load_vel,
                                                            create_hdf5):
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

    readFunc = read_structured_velocity_hdf5(create_hdf5, addValTop=(0, 0, 0),
                                             excludeBot=nPointsY, interpValBot=1)

    [uXR, uYR, uZR] = readFunc(0)

    assert np.all(uX == uXR)
    assert np.all(uY == uYR)
    assert np.all(uZ == uZR)

