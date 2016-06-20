import eddylicious
from eddylicious.writers.hdf5_writers import *
import numpy as np
import h5py
import pytest
from os import path


def test_point_writer(tmpdir):
    dsvDir = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180", "dsv_output")

    pointsY = np.load(path.join(dsvDir, "pointsY.npy"))
    pointsZ = np.load(path.join(dsvDir, "pointsZ.npy"))
    xVal = 0.0
    writeFile = tmpdir.join("test.hdf5")
    writePath = writeFile.strpath
    dbFile = h5py.File(writePath, 'a')
    write_points_to_hdf5(dbFile, pointsY, pointsZ, xVal)

    writtenPoints = dbFile["points"]

    assert np.all(writtenPoints[:, 0] ==
                  xVal*np.ones(writtenPoints[:, 0].shape))
    assert np.all(writtenPoints[:, 1] ==
                  pointsY.reshape((pointsY.size, -1), order='F')[:, 0])

    assert np.all(writtenPoints[:, 2] ==
                  pointsZ.reshape((pointsY.size, -1), order='F')[:, 0])


def test_velocity_writer(tmpdir):
    dsvDir = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180", "dsv_output")

    uX = np.load(path.join(dsvDir, "1000.05", "uX.npy"))
    uY = np.load(path.join(dsvDir, "1000.05", "uY.npy"))
    uZ = np.load(path.join(dsvDir, "1000.05", "uZ.npy"))

    size = 3
    iteration = 1

    writeFile = tmpdir.join("test.hdf5")
    writePath = writeFile.strpath
    dbFile = h5py.File(writePath, 'a')
    dbFile.create_dataset("time", data=np.ones((size, 1)))
    dbFile.create_dataset("velocity", (size, uX.size, 3), dtype=np.float)
    write_velocity_to_hdf5(dbFile, 0.1, uX, uY, uZ, iteration)

    uX = uX.reshape((uX.size, -1), order='F')[:, 0]
    uY = uY.reshape((uY.size, -1), order='F')[:, 0]
    uZ = uZ.reshape((uZ.size, -1), order='F')[:, 0]

    assert np.all(dbFile["velocity"][iteration, :, 0] == uX)
    assert np.all(dbFile["velocity"][iteration, :, 1] == uY)
    assert np.all(dbFile["velocity"][iteration, :, 2] == uZ)


def test_velocity_writer_iter_larger_than_total_size(tmpdir):
    dsvDir = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                       "channel_flow_180", "dsv_output")

    uX = np.load(path.join(dsvDir, "1000.05", "uX.npy"))
    uY = np.load(path.join(dsvDir, "1000.05", "uY.npy"))
    uZ = np.load(path.join(dsvDir, "1000.05", "uZ.npy"))

    size = 3
    iteration = 3

    writeFile = tmpdir.join("test.hdf5")
    writePath = writeFile.strpath
    dbFile = h5py.File(writePath, 'a')
    dbFile.create_dataset("time", data=np.ones((size, 1)))
    dbFile.create_dataset("velocity", (size, uX.size, 3), dtype=np.float)

    with pytest.raises(ValueError):
        write_velocity_to_hdf5(dbFile, 0.1, uX, uY, uZ, iteration)
