import os
import numpy as np
import h5py as h5py
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_u_from_foamfile
import argparse


# Define the command-line arguments
parser = argparse.ArgumentParser(
            description="A script for converting a database stored \
                    as a collection of foamFile-formatted files to \
                    a single hdf5 file. \
                    Produces the following files: \
                    uMean, uPrimeX, uPrimeY, uPrimeZ, pointsY, pointsZ."
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

# Get the mean profile
uMean = np.append(np.zeros((1, 1)), np.genfromtxt(uMeanFile)[:, 1])

nPointsY = uMean.size

# Read in the points
[pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
    os.path.join(dataDir, times[0], surfaceName, "faceCentres"),
    nPointsY=nPointsY)

[nPointsY, nPointsZ] = pointsY.shape

# Allocate arrays for the fluctuations
if os.path.isfile(fileName):
    print "HDF5 file already exsists. It it will be overwritten."
    os.remove(fileName)

dbFile = h5py.File(fileName, 'a')

pointsGroup = dbFile.create_group("points")
velocityGroup = dbFile.create_group("velocity")

pointsGroup.create_dataset("pointsY", data=pointsY)
pointsGroup.create_dataset("pointsZ", data=pointsZ)

velocityGroup.create_dataset("uMean", data=uMean)


uPrimeX = velocityGroup.create_dataset("uPrimeX", (len(times),
                                                   pointsY.shape[0],
                                                   pointsY.shape[1]))
uPrimeY = velocityGroup.create_dataset("uPrimeY", (len(times),
                                                   pointsY.shape[0],
                                                   pointsY.shape[1]))
uPrimeZ = velocityGroup.create_dataset("uPrimeZ", (len(times),
                                                   pointsY.shape[0],
                                                   pointsY.shape[1]))

dbFile.attrs["nPointsY"] = pointsY.shape[0]
dbFile.attrs["nPointsZ"] = pointsY.shape[1]
dbFile.attrs["nPoints"] = pointsY.size

printCounter = 0

# Read in the fluctuations
for timeI in xrange(len(times)):
    printCounter += 1

    if printCounter == 10:
        printCounter = 0
        print "Read in", timeI, "time-iterations out of", len(times), "."

    # Read in U
    [uX, uY, uZ] = read_u_from_foamfile(os.path.join(dataDir,
                                                     times[timeI],
                                                     surfaceName,
                                                     "vectorField", "U"),
                                        nPointsY, nPointsZ, yInd, zInd)

    uPrimeX[timeI, :, :] = uX - uMean[:, np.newaxis]
    uPrimeY[timeI, :, :] = uY
    uPrimeZ[timeI, :, :] = uZ

print pointsGroup["pointsY"]
dbFile.close()
