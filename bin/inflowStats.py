from __future__ import print_function
from __future__ import division
import numpy as np
import os as os
import h5py
import argparse
from mpi4py import MPI
from eddylicious.generators.helper_functions import chunks_and_offsets


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()


# Define the command-line arguments
    parser = argparse.ArgumentParser(
                description="A utility for calculating statistics \
                        of a inflow field database stored in hdf5 format.")

    parser.add_argument('--database',
                        type=str,
                        help='The HDF5 file with the database.',
                        required=True)
    parser.add_argument('--writePath',
                        type=str,
                        help='The location where to write the \
                            computed results.',
                        required=True)

    args = parser.parse_args()

    readPath = args.readPath
    writeDir = args.writePath

# Open the hdf5 database
    if rank == 0:
        print("Opening the database")

    dbFile = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
    times = dbFile['time'][()]
    points = dbFile['points'][()]
    points = points[:, 1:]

    size = len(times)

    uMean = np.zeros((points.shape[0], 3))
    uSquaredMean = np.zeros((points.shape[0], 3))

    if rank == 0:
        print("Calculating the statistics")

    [chunks, offsets] = chunks_and_offsets(nProcs, size)

    for i in range(chunks[rank]):
        if rank == 0 and (np.mod(i, int(chunks[rank]/10)) == 0):

            print("Computed about " + str(int(i/chunks[rank]*100)) + "%")

        position = offsets[rank] + i

        uMean += dbFile['velocity'][position, :, :]
        uSquaredMean += dbFile['velocity'][position, :, :]**2

    comm.Barrier()
    dbFile.close()
    if rank == 0:
        print("Done")

    uMean = comm.gather(uMean, root=0)
    uSquaredMean = comm.gather(uSquaredMean, root=0)

    if rank == 0:
        for i in range(nProcs-1):
            uMean[0] += uMean[i+1]
            uSquaredMean[0] += uSquaredMean[i+1]

        uMean = uMean[0]/size
        uSquaredMean = uSquaredMean[0]/size

        uPrime2Mean = uSquaredMean - uMean**2

        print("Reshaping and averaging")
        # Sort along y first
        yInd = np.argsort(points[:, 0])
        points[:, 0] = points[yInd, 0]
        points[:, 1] = points[yInd, 1]
        uMean[:, 0] = uMean[yInd, 0]
        uMean[:, 1] = uMean[yInd, 1]
        uMean[:, 2] = uMean[yInd, 2]
        uPrime2Mean[:, 0] = uPrime2Mean[yInd, 0]
        uPrime2Mean[:, 1] = uPrime2Mean[yInd, 1]
        uPrime2Mean[:, 2] = uPrime2Mean[yInd, 2]

        # Find the number of points along z
        nPointsZ = 0
        for i in range(points[:, 0].size):
            if points[i, 0] == points[0, 0]:
                nPointsZ += 1
            else:
                break

        # Reshape into a 2d array
        pointsY = np.copy(np.reshape(points[:, 0], (-1, nPointsZ)))
        pointsZ = np.copy(np.reshape(points[:, 1], (-1, nPointsZ)))
        uMeanX = np.copy(np.reshape(uMean[:, 0], (-1, nPointsZ)))
        uMeanY = np.copy(np.reshape(uMean[:, 1], (-1, nPointsZ)))
        uMeanZ = np.copy(np.reshape(uMean[:, 2], (-1, nPointsZ)))
        uPrime2MeanXX = np.copy(np.reshape(uPrime2Mean[:, 0],
                                           (-1, nPointsZ)))
        uPrime2MeanYY = np.copy(np.reshape(uPrime2Mean[:, 1],
                                           (-1, nPointsZ)))
        uPrime2MeanZZ = np.copy(np.reshape(uPrime2Mean[:, 2],
                                           (-1, nPointsZ)))

        # For each y order the points in z

        zInd = np.zeros(pointsZ.shape, dtype=np.int)

        for i in range(pointsZ.shape[0]):
            zInd[i, :] = np.argsort(pointsZ[i, :])
            pointsZ[i, :] = pointsZ[i, zInd[i, :]]
            uMeanX[i, :] = uMeanX[i, zInd[i, :]]
            uMeanY[i, :] = uMeanY[i, zInd[i, :]]
            uMeanZ[i, :] = uMeanZ[i, zInd[i, :]]
            uPrime2MeanXX[i, :] = uPrime2MeanXX[i, zInd[i, :]]
            uPrime2MeanYY[i, :] = uPrime2MeanYY[i, zInd[i, :]]
            uPrime2MeanZZ[i, :] = uPrime2MeanZZ[i, zInd[i, :]]

        y = pointsY[:, 0]

        # Average along Z
        uMeanX = np.mean(uMeanX, axis=1)
        uMeanY = np.mean(uMeanY, axis=1)
        uMeanZ = np.mean(uMeanZ, axis=1)

        uPrime2MeanXX = np.mean(uPrime2MeanXX, axis=1)
        uPrime2MeanYY = np.mean(uPrime2MeanYY, axis=1)
        uPrime2MeanZZ = np.mean(uPrime2MeanZZ, axis=1)

        print("Outputting data")

        if not os.path.exists(writeDir):
            os.makedirs(writeDir)

        np.savetxt(os.path.join(writeDir, "y"), y)
        np.savetxt(os.path.join(writeDir, "uMeanX"), uMeanX)
        np.savetxt(os.path.join(writeDir, "uMeanY"), uMeanY)
        np.savetxt(os.path.join(writeDir, "uMeanZ"), uMeanZ)
        np.savetxt(os.path.join(writeDir, "uPrime2MeanXX"), uPrime2MeanXX)
        np.savetxt(os.path.join(writeDir, "uPrime2MeanYY"), uPrime2MeanYY)
        np.savetxt(os.path.join(writeDir, "uPrime2MeanZZ"), uPrime2MeanZZ)
        np.savetxt(os.path.join(writeDir, "y"), y)

if __name__ == "__main__":
    main()

