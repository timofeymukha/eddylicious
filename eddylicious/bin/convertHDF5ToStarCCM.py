from __future__ import print_function
from __future__ import division
import os
import numpy as np
import h5py
from mpi4py import MPI
from eddylicious.generators.helper_functions import chunks_and_offsets
import argparse


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()


# Define the command-line arguments
    parser = argparse.ArgumentParser(
                description="A utility for converting a database stored \
                            as a HDF5 file into a set of CSV files \
                            supported by STAR-CCM+."
                                    )

    parser.add_argument('--database',
                        type=str,
                        help='The hdf5 to convert.',
                        required=True)

    parser.add_argument('--writepath',
                        type=str,
                        help='Directory to which the files will be written.',
                        required=True)

    args = parser.parse_args()

    database = args.database
    writePath = args.writepath

# Open the hdf5 file
    dbFile = h5py.File(database, 'r', driver='mpio', comm=MPI.COMM_WORLD)

# Grab the existing times
    times = dbFile["velocity"]["times"]

# Get the mean profile and append zeros
    uMeanX = dbFile["velocity"]["uMeanX"]
    uMeanY = dbFile["velocity"]["uMeanY"]

# Read in the points
    pointsY = dbFile["pointsY"]["pointsY"]
    pointsZ = dbFile["pointsZ"]["pointsZ"]
    np.reshape(pointsY, (pointsY.size, -1), order='F')
    np.reshape(pointsZ, (pointsZ.size, -1), order='F')

    [nPointsY, nPointsZ] = pointsY.shape

    assert nPointsY == uMean.shape[0]

    fileHeader = "X Y Z uX uY uZ"

# Distribute time frame between processors
    [chunks, offsets] = chunks_and_offsets(nProcs, len(times))

    for i in range(chunks[rank]):
        if rank == 0 and (np.mod(i, int(chunks[rank]/20)) == 0):
            print("Converted about " + str(i/chunks[rank]*100)+"%")

        position = offsets[rank] + i
        # Read in U
        uX = dbFile["velocity"]["uX"][position, :, :]
        uY = dbFile["velocity"]["uY"][position, :, :]
        uZ = dbFile["velocity"]["uZ"][position, :, :]
        np.reshape(uY, (uY.size, -1), order='F')
        np.reshape(uX, (uX.size, -1), order='F')
        np.reshape(uZ, (uY.size, -1), order='F')

        np.savetxt(str(position)+".csv",
                   np.concatenate((np.zeros(pointsY.shape), pointsY, pointsZ,
                                   uX, uY, uZ)),
                   header=fileHeader)

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()
    dbFile.close()
    if rank == 0:
        print("Done")

if __name__ == " __main__":
    main()
