import numpy as np
import os as os
import matplotlib.pyplot as plt
import h5py as h5py
import argparse


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

print "Opening the database"
dbFile = h5py.File(readPath, 'r')

pointsY = dbFile["points"]["pointsY"]
pointsZ = dbFile["points"]["pointsZ"]
uX = dbFile["velocity"]["uX"]
uY = dbFile["velocity"]["uY"]
uZ = dbFile["velocity"]["uZ"]

size = len(uX[:,0])
nPointsY = pointsY.shape[0]
nPointsZ = pointsY.shape[1]

uMean = np.zeros((nPointsY, nPointsZ, 3))
uSquaredMean = np.zeros((nPointsY, nPointsZ, 3))

print "Calculating the statistics"

for i in xrange(size):
    uMean[:, :, 0] += uX[i, :, :]
    uMean[:, :, 1] += uY[i, :, :]
    uMean[:, :, 2] += uZ[i, :, :]
    uSquaredMean[:, :, 0] += uX[i, :, :]**2
    uSquaredMean[:, :, 1] += uY[i, :, :]**2
    uSquaredMean[:, :, 2] += uZ[i, :, :]**2


uMean /= size
uSquaredMean /= size

uPrime2Mean = uSquaredMean - uMean**2

# Average along Z
uMean = np.mean(uMean, axis=1)

uPrime2Mean = np.mean(uPrime2Mean, axis=1)

print "Outputting figures and data"

if not os.path.exists(writeDir):
    os.makedirs(writeDir)

np.savetxt(os.path.join(writeDir, "y"), pointsY[:, 0])
np.savetxt(os.path.join(writeDir, "uMeanX"), uMean[:,0])
np.savetxt(os.path.join(writeDir, "uMeanY"), uMeanY[:,1])
np.savetxt(os.path.join(writeDir, "uMeanZ"), uMeanZ[:,2])
np.savetxt(os.path.join(writeDir, "uPrime2MeanX"), uPrime2Mean[:,0])
np.savetxt(os.path.join(writeDir, "uPrime2MeanY"), uPrime2Mean[:,1])
np.savetxt(os.path.join(writeDir, "uPrime2MeanZ"), uPrime2Mean[:,2])
np.savetxt(os.path.join(writeDir, "y"), pointsY[:, 0])
