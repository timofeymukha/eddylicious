import numpy as np
import os as os
import matplotlib.pyplot as plt
import h5py as h5py
import argparse
from mpi4py import MPI
from eddylicious.generators.helper_functions import chunks_and_offsets

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nProcs = comm.Get_size()

# Define the command-line arguments
parser = argparse.ArgumentParser(
            description="A script for plotting statistics \
                    of a precursor database.")

parser.add_argument('--readPath',
                    type=str,
                    help='The HDF5 file with the databse.',
                    required=True)
parser.add_argument('--writePath',
                    type=str,
                    help='The location where to write the \
                          produced plots.',
                    required=True)

args = parser.parse_args()

readPath = args.readPath
writeDir = args.writePath


#  Open the hdf5 database

if rank == 0:
    print "Opening the database"
dbFile = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)

pointsY = dbFile["points"]["pointsY"]
pointsZ = dbFile["points"]["pointsZ"]

size = len(uX[:,0])
nPointsY = pointsY.shape[0]
nPointsZ = pointsY.shape[1]

uMean = np.zeros((nPointsY, nPointsZ, 3))
uSquaredMean = np.zeros((nPointsY, nPointsZ, 3))

if rank == 0:
    print "Calculating the statistics"

[chunks, offsets] = chunks_and_offsets(nProcs, size)

for i in xrange(chunks[rank]):
    if rank == 0:
        print "Computed about", i/float(chunks[rank])*100, "%"

    position = offsets[rank] + i

    uMean[:, :, 0] += dbFile["velocity"]["uX"][position, :, :]
    uMean[:, :, 1] += dbFile["velocity"]["uY"][position, :, :]
    uMean[:, :, 2] += dbFile["velocity"]["uZ"][position, :, :]
    uSquaredMean[:, :, 0] += dbFile["velocity"]["uX"][position, :, :]**2
    uSquaredMean[:, :, 1] += dbFile["velocity"]["uY"][position, :, :]**2
    uSquaredMean[:, :, 2] += dbFile["velocity"]["uZ"][position, :, :]**2

uMean = comm.gather(uMean, root=0)
uSquaredMean = comm.gather(uSquaredMean, root=0)

if rank == 0:
    for i in xrange(nProcs-1):
        uMean[0] += uMean[i+1]
        uSquaredMean[0] += uSquaredMean[i+1]

    uMean = uMean[0]/size
    uSquaredMean = uSquaredMean[0]/size

    uPrime2Mean = uSquaredMean - uMean**2

# Average along Z
    uMean = np.mean(uMean, axis=1)
    uPrime2Mean = np.mean(uPrime2Mean, axis=1)

    print "Outputting data"

    if not os.path.exists(writeDir):
        os.makedirs(writeDir)

    np.savetxt(os.path.join(writeDir, "y"), pointsY[:, 0])
    np.savetxt(os.path.join(writeDir, "uMeanX"), uMean[:,0])
    np.savetxt(os.path.join(writeDir, "uMeanY"), uMean[:,1])
    np.savetxt(os.path.join(writeDir, "uMeanZ"), uMean[:,2])
    np.savetxt(os.path.join(writeDir, "uPrime2MeanX"), uPrime2Mean[:,0])
    np.savetxt(os.path.join(writeDir, "uPrime2MeanY"), uPrime2Mean[:,1])
    np.savetxt(os.path.join(writeDir, "uPrime2MeanZ"), uPrime2Mean[:,2])
    np.savetxt(os.path.join(writeDir, "y"), pointsY[:, 0])
