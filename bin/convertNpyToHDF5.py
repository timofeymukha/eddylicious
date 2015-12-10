import os as os
import numpy as np
import h5py as h5py


dbDir = "dbtest"

# Load database in npy format
print "Loading database.."
uMean = np.load(os.path.join(dbDir, "uMean.npy"))
pointsY = np.load(os.path.join(dbDir, "pointsY.npy"))
pointsZ = np.load(os.path.join(dbDir, "pointsZ.npy"))

uPrimeX = np.load(os.path.join(dbDir, "uPrimeX.npy"))
uPrimeY = np.load(os.path.join(dbDir, "uPrimeY.npy"))
uPrimeZ = np.load(os.path.join(dbDir, "uPrimeZ.npy"))

# Create the hdf5 database
dbFile = h5py.File('testDb.hdf5', 'a')

pointsGroup = dbFile.create_group("points")
velocityGroup = dbFile.create_group("velocity")

pointsGroup.create_dataset("pointsY", data=pointsY)
pointsGroup.create_dataset("pointsZ", data=pointsY)

velocityGroup.create_dataset("uMean", data=uMean)
velocityGroup.create_dataset("uPrimeX", data=uPrimeX, compression="gzip")
velocityGroup.create_dataset("uPrimeY", data=uPrimeY, compression="gzip")
velocityGroup.create_dataset("uPrimeZ", data=uPrimeZ, compression="gzip")

# Add attributes
dbFile.attrs["nu"] = 0.001973
dbFile.attrs["uTau"] = 1.0
dbFile.attrs["delta"] = 1.0
dbFile.attrs["deltaT"] = 0.0008

pointsGroup.attrs["nPointsY"] = pointsY.shape[0]
pointsGroup.attrs["nPointsZ"] = pointsY.shape[1]
pointsGroup.attrs["nPoints"] = pointsY.size
