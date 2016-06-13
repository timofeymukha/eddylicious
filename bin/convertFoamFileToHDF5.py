from __future__ import print_function
from __future__ import division
import os
import numpy as np
import h5py
from mpi4py import MPI
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_from_foamfile
from eddylicious.generators.helper_functions import chunks_and_offsets
import argparse


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()


# Define the command-line arguments
    parser = argparse.ArgumentParser(
                description="A utility for converting a database stored \
                            as a collection of foamFile-formatted files to \
                            a single HDF5 file."
                                    )

    parser.add_argument('--precursor',
                        type=str,
                        help='The location of the precusor case.',
                        required=True)
    parser.add_argument('--surface',
                        type=str,
                        help='The name of the surface that contains the data.',
                        required=True)
    parser.add_argument('--filename',
                        type=str,
                        help='The name hdf5 file to create.',
                        required=True)

    parser.add_argument('--umean',
                        type=str,
                        help='The file containing the mean velocity profile.',
                        required=True)

    args = parser.parse_args()

    precursorCaseDir = args.precursor
    surfaceName = args.surface
    uMeanFile = args.umean
    fileName = args.filename

    dataDir = os.path.join(precursorCaseDir, "postProcessing",
                           "sampledSurface")

# Grab the existing times and sort
    times = os.listdir(dataDir)
    times = np.sort(times)

# Get the mean profile and append zeros
    uMean = np.genfromtxt(uMeanFile)
    uMeanX = uMean[:, 1]
    if uMean.shape[1] == 3:
        uMeanY = uMean[:, 2]
    else:
        uMeanY = np.zeros(uMeanX.shape)

    y = uMean[:, 0]


# Read in the points
    [pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
        os.path.join(dataDir, times[0], surfaceName, "faceCentres"),
        addValBot=y[0], addValTop=y[-1])

    [nPointsY, nPointsZ] = pointsY.shape

    assert nPointsY == uMean.shape[0]


# Allocate arrays for the fluctuations
    if rank == 0:
        if os.path.isfile(fileName):
            print("HDF5 file already exists. It it will be overwritten.")
            os.remove(fileName)

    dbFile = h5py.File(fileName, 'a', driver='mpio', comm=MPI.COMM_WORLD)

    pointsGroup = dbFile.create_group("points")
    velocityGroup = dbFile.create_group("velocity")

    pointsGroup.create_dataset("pointsY", data=pointsY)
    pointsGroup.create_dataset("pointsZ", data=pointsZ)

    velocityGroup.create_dataset("uMeanX", data=uMeanX)
    velocityGroup.create_dataset("uMeanY", data=uMeanY)

    velocityGroup.create_dataset("times", data=[float(times[i])
                                                for i in range(times.size)])

    uX = velocityGroup.create_dataset("uX", (len(times),
                                             pointsY.shape[0],
                                             pointsY.shape[1]),
                                      dtype=np.float64)
    uY = velocityGroup.create_dataset("uY", (len(times),
                                             pointsY.shape[0],
                                             pointsY.shape[1]),
                                      dtype=np.float64)
    uZ = velocityGroup.create_dataset("uZ", (len(times),
                                             pointsY.shape[0],
                                             pointsY.shape[1]),
                                      dtype=np.float64)

    dbFile.attrs["nPointsY"] = pointsY.shape[0]
    dbFile.attrs["nPointsZ"] = pointsY.shape[1]
    dbFile.attrs["nPoints"] = pointsY.size

    [chunks, offsets] = chunks_and_offsets(nProcs, len(times))

    readFunc = read_velocity_from_foamfile(dataDir, surfaceName,
                                           nPointsZ, yInd, zInd,
                                           addValBot=(uMeanX[0], uMeanY[0], 0),
                                           addValTop=(uMeanX[-1], uMeanY[-1],
                                                      0))

# Read in the fluctuations
    for i in range(chunks[rank]):
        if rank == 0 and (np.mod(i, int(chunks[rank]/20)) == 0):
            print("Converted about " + str(i/chunks[rank]*100)+"%")

        position = offsets[rank] + i
        # Read in U
        [uXVal, uYVal, uZVal] = readFunc(times[position])

        uX[position, :, :] = uXVal
        uY[position, :, :] = uYVal
        uZ[position, :, :] = uZVal

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()
    dbFile.close()
    if rank == 0:
        print("Done")

if __name__ == " __main__":
    main()
