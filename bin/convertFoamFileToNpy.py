import os
import numpy as np
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_u_from_foamfile
import argparse


# Define the command-line arguments
parser = argparse.ArgumentParser(
            description="A script for converting a database stored \
                    as a collection of foamFile-formatted files to \
                    npy files. \
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
parser.add_argument('--writePath',
                    type=str,
                    help='The location where to write the \
                          produced numpy arrays.',
                    required=True)
parser.add_argument('--uMeanFile',
                    type=str,
                    help='The file containg the mean velocity profile. \
                          The file is assumed to have to columns, \
                          one with y coordinates, and the other one \
                          with the values of mean velocity.',
                    required=True)

args = parser.parse_args()

precursorCaseDir = args.precursorPath
surfaceName = args.surfaceName
writeDir = args.writePath
uMeanFile = args.uMeanFile


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
uPrimeX = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeY = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeZ = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))


printCounter = 0

# Read in the fluctuations
for timeI in xrange(len(times)):
    printCounter += 1

    if (printCounter == 100):
        printCounter = 0
        print "Read in", timeI, "time-iterations out of", len(times), "."

    # Read in U
    [uX, uY, uZ] = read_u_from_foamfile(
        os.path.join(dataDir, times[timeI], surfaceName, "vectorField", "U"),
        nPointsY, nPointsZ,
        yInd,
        zInd)

    uPrimeX[:, :, timeI] = uX - uMean[:, np.newaxis]
    uPrimeY[:, :, timeI] = uY
    uPrimeZ[:, :, timeI] = uZ

uPrimeX = uPrimeX[:, :, 1:]
uPrimeY = uPrimeY[:, :, 1:]
uPrimeZ = uPrimeZ[:, :, 1:]

# Write to disk ad npy files
np.save(os.path.join(writeDir, "uPrimeX"), uPrimeX)
np.save(os.path.join(writeDir, "uPrimeY"), uPrimeY)
np.save(os.path.join(writeDir, "uPrimeZ"), uPrimeZ)
np.save(os.path.join(writeDir, "pointsY"), pointsY)
np.save(os.path.join(writeDir, "pointsZ"), pointsZ)
np.save(os.path.join(writeDir, "uMean"), uMean)
