import os
import numpy as np
import h5py as h5py
from mpi4py import MPI
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_u_from_foamfile
from eddylicious.generators.helper_functions import chunks_and_offsets
import argparse


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nProcs = comm.Get_size()


# Define the command-line arguments
parser = argparse.ArgumentParser(
            description="A script for converting a database stored \
                        as a collection of foamFile-formatted files to \
                        a single hdf5 file. \
                        The file contains two groups: points and velocity. \
                        The points group contains the datasets pointsY and \
                        pointsZ. \
                        The velocity group contains the datasets: \
                        uMean, uX, uY and uZ. \
                        Also, the attributes nPointsY, nPointsZ and nPoints \
                        are added to the root of the file."
                                )

parser.add_argument('--precursorPath',
                    type=str,
                    help='The location of the precusor case',
                    required=True)
parser.add_argument('--surfaceName',
                    type=str,
                    help='The name of the surface that contains the data.',
                    required=True)
parser.add_argument('--fileName',
                    type=str,
                    help='The location where to write the \
                          produced hdf5 file.',
                    required=True)
parser.add_argument('--uMeanFile',
                    type=str,
                    help='The file containing the mean velocity profile. \
                          The file is assumed to have two columns, \
                          one with y coordinates, and the other one \
                          with the values of mean velocity.',
                    required=True)

args = parser.parse_args()

precursorCaseDir = args.precursorPath
surfaceName = args.surfaceName
uMeanFile = args.uMeanFile
fileName = args.fileName


dataDir = os.path.join(precursorCaseDir, "postProcessing", "sampledSurface")

# Grab the existing times and sort
times = os.listdir(dataDir)
times = np.sort(times)

# Get the mean profile and append zeros
uMean = np.append(np.zeros((1, 1)), np.genfromtxt(uMeanFile)[:, 1])
uMean = np.append(uMean, np.zeros((1, 1)))

# Read in the points
[pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
    os.path.join(dataDir, times[0], surfaceName, "faceCentres"),
    addValBot=0, addValTop=2)

[nPointsY, nPointsZ] = pointsY.shape

assert nPointsY == uMean.size

# Allocate arrays for the fluctuations
if rank == 0:
    if os.path.isfile(fileName):
        print "HDF5 file already exsists. It it will be overwritten."
        os.remove(fileName)

dbFile = h5py.File(fileName, 'a', driver='mpio', comm=MPI.COMM_WORLD)

pointsGroup = dbFile.create_group("points")
velocityGroup = dbFile.create_group("velocity")

pointsGroup.create_dataset("pointsY", data=pointsY)
pointsGroup.create_dataset("pointsZ", data=pointsZ)

velocityGroup.create_dataset("uMean", data=uMean)

velocityGroup.create_dataset("times", data=[float(times[i])
                                            for i in xrange(times.size)])

uX = velocityGroup.create_dataset("uX", (len(times),
                                         pointsY.shape[0],
                                         pointsY.shape[1]))
uY = velocityGroup.create_dataset("uY", (len(times),
                                         pointsY.shape[0],
                                         pointsY.shape[1]))
uZ = velocityGroup.create_dataset("uZ", (len(times),
                                         pointsY.shape[0],
                                         pointsY.shape[1]))

dbFile.attrs["nPointsY"] = pointsY.shape[0]
dbFile.attrs["nPointsZ"] = pointsY.shape[1]
dbFile.attrs["nPoints"] = pointsY.size

[chunks, offsets] = chunks_and_offsets(nProcs, len(times))

# Read in the fluctuations
for i in xrange(chunks[rank]):
    if rank == 0:
        print "Converted about", i/float(chunks[rank])*100, "%"

    position = offsets[rank] + i
    # Read in U
    [uXVal, uYVal, uZVal] = read_u_from_foamfile(os.path.join(dataDir,
                                                              times[position],
                                                              surfaceName,
                                                              "vectorField",
                                                              "U"), nPointsY,
                                                 nPointsZ, yInd, zInd,
                                                 addValBot=0,
                                                 addValTop=0,
                                                 interpolate=False)

    uX[position, :, :] = uXVal
    uY[position, :, :] = uYVal
    uZ[position, :, :] = uZVal

if rank == 0:
    print "Wait for it.. Other processes finishing up."
