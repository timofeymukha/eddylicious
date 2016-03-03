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

times = dbFile['time'][()]
points = dbFile['points'][()]
points = points[:, 1:]
velocity = dbFile['velocity']

size = len(times)

uMean = np.zeros((points.shape[0], 3))
uSquaredMean = np.zeros((points.shape[0], 3))

print "Calculating the statistics"

for i in xrange(size):
    uMean += velocity[i, :, :]
    uSquaredMean += velocity[i, :, :]**2


uMean /= size
uSquaredMean /= size

uPrime2Mean = uSquaredMean - uMean**2

print "Reshaping and averaging"
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
for i in xrange(points[:, 0].size):
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
uPrime2MeanX = np.copy(np.reshape(uPrime2Mean[:, 0],
                                  (-1, nPointsZ)))
uPrime2MeanY = np.copy(np.reshape(uPrime2Mean[:, 1],
                                  (-1, nPointsZ)))
uPrime2MeanZ = np.copy(np.reshape(uPrime2Mean[:, 2],
                                  (-1, nPointsZ)))

# For each y order the points in z

zInd = np.zeros(pointsZ.shape, dtype=np.int)

for i in xrange(pointsZ.shape[0]):
    zInd[i, :] = np.argsort(pointsZ[i, :])
    pointsZ[i, :] = pointsZ[i, zInd[i, :]]
    uMeanX[i, :] = uMeanX[i, zInd[i, :]]
    uMeanY[i, :] = uMeanY[i, zInd[i, :]]
    uMeanZ[i, :] = uMeanZ[i, zInd[i, :]]
    uPrime2MeanX[i, :] = uPrime2MeanX[i, zInd[i, :]]
    uPrime2MeanY[i, :] = uPrime2MeanY[i, zInd[i, :]]
    uPrime2MeanZ[i, :] = uPrime2MeanZ[i, zInd[i, :]]

y = pointsY[:, 0]

# Average along Z
uMeanX = np.mean(uMeanX, axis=1)
uMeanY = np.mean(uMeanY, axis=1)
uMeanZ = np.mean(uMeanZ, axis=1)

uPrime2MeanX = np.mean(uPrime2MeanX, axis=1)
uPrime2MeanY = np.mean(uPrime2MeanY, axis=1)
uPrime2MeanZ = np.mean(uPrime2MeanZ, axis=1)

print "Outputting figures and data"

if not os.path.exists(writeDir):
    os.makedirs(writeDir)

np.savetxt(os.path.join(writeDir, "y"), y)
np.savetxt(os.path.join(writeDir, "uMeanX"), uMeanX)
np.savetxt(os.path.join(writeDir, "uMeanY"), uMeanY)
np.savetxt(os.path.join(writeDir, "uMeanZ"), uMeanZ)
np.savetxt(os.path.join(writeDir, "uPrime2MeanX"), uPrime2MeanX)
np.savetxt(os.path.join(writeDir, "uPrime2MeanY"), uPrime2MeanY)
np.savetxt(os.path.join(writeDir, "uPrime2MeanZ"), uPrime2MeanZ)
np.savetxt(os.path.join(writeDir, "y"), y)

#plt.figure()
#plt.plot(y, uMeanX)
#plt.xlabel("y")
#plt.ylabel(r"$<u_x>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uMeanX.png"), format="png")

#plt.figure()
#plt.plot(y, uMeanY)
#plt.xlabel("y")
#plt.ylabel(r"$<u_y>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uMeanY.png"), format="png")

#plt.figure()
#plt.plot(y, uMeanZ)
#plt.xlabel("y")
#plt.ylabel(r"$<u_z>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uMeanZ.png"), format="png")

#plt.figure()
#plt.plot(y, uPrime2MeanX)
#plt.xlabel("y")
#plt.ylabel(r"$<u'_xu'_x>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uPrime2X.png"), format="png")

#plt.figure()
#plt.plot(y, uPrime2MeanY)
#plt.xlabel("y")
#plt.ylabel(r"$<u'_yu'_y>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uPrime2Y.png"), format="png")

#plt.figure()
#plt.plot(y, uPrime2MeanZ)
#plt.xlabel("y")
#plt.ylabel(r"$<u'_zu'_z>$")
#plt.grid()
#plt.savefig(os.path.join(writeDir, "uPrime2Z.png"), format="png")


