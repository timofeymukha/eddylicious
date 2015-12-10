import os
import numpy as np
from ..readers.foamfile_readers import read_points_from_foamfile
from ..readers.foamfile_readers import read_u_from_foamfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--precursorPath',
                    type=str,
                    help='The location of the precusor case')
parser.add_argument('--writePath',
                    type=str,
                    help='The location where to write the produced numpy arrays')

args = parser.parse_args()


precursorCaseDir = "../channelPronkRefined"

surfaceName = "inletSurface"

dataDir = os.path.join(precursorCaseDir, "postProcessing", "sampledSurface")


# Grab the existing times and sort
times = os.listdir(dataDir)

# Get the mean profile
uMean = np.append(np.zeros((1, 1)), np.genfromtxt(
                os.path.join(precursorlCaseDir,
                             "postProcessing",
                             "collapsedFields",
                             "240", "UMean_X.xy"))[:, 1])

nPointsY = uMean.size

# Read grid for the recycling plane
[pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
    os.path.join(dataDir, times[0], surfaceName, "faceCentres"),
    nPointsY=nPointsY)

[nPointsY, nPointsZ] = pointsY.shape

uPrimeX = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeY = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeZ = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))

times = np.sort(times)

for timeI in xrange(len(times)):

    # Read U data
    [uX, uY, uZ] = read_u_from_foamfile(
        os.path.join(dataDir, times[timeI], surfaceName, "vectorField", "U"),
        nPointsY,
        yInd,
        zInd)

    uPrimeX[:, :, timeI] = uX - uMean[:, np.newaxis]
    uPrimeY[:, :, timeI] = uY
    uPrimeZ[:, :, timeI] = uZ

    print timeI

uPrimeX = uPrimeX[:, :, 1:]
uPrimeY = uPrimeY[:, :, 1:]
uPrimeZ = uPrimeZ[:, :, 1:]

np.save(os.path.join(writeDir, "uPrimeX"), uPrimeX)
np.save(os.path.join(writeDir, "uPrimeY"), uPrimeY)
np.save(os.path.join(writeDir, "uPrimeZ"), uPrimeZ)
np.save(os.path.join(writeDir, "pointsY"), pointsY)
np.save(os.path.join(writeDir, "pointsZ"), pointsZ)
np.save("dbtest/uMean", uMean)
