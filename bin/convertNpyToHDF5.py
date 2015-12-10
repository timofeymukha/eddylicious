import os as os
import numpy as np
import h5py as h5py
import argparse


# Define the command-line arguments
parser = argparse.ArgumentParser(
            description="A script for converting a database stored \
                    as a collection of npy files to a HDF5 file.")

parser.add_argument('--readPath',
                    type=str,
                    help='The directory holding the npuy files,',
                    required=True)
parser.add_argument('--writePath',
                    type=str,
                    help='The location where to write the \
                          produced HDF5 file.',
                    required=True)
parser.add_argument('--fileName',
                    type=str,
                    help='The name of the produced HDF5 file.',
                    required=True)

args = parser.parse_args()

readPath = args.readPath
writeDir = args.writePath


# Load database in npy format
print "Loading database.."
uMean = np.load(os.path.join(readPath, "uMean.npy"))
pointsY = np.load(os.path.join(readPath, "pointsY.npy"))
pointsZ = np.load(os.path.join(readPath, "pointsZ.npy"))

uPrimeX = np.load(os.path.join(readPath, "uPrimeX.npy"))
uPrimeY = np.load(os.path.join(readPath, "uPrimeY.npy"))
uPrimeZ = np.load(os.path.join(readPath, "uPrimeZ.npy"))

print "Done. Creating the HDF5 file."

# Create the hdf5 database
dbFile = h5py.File(args.fileName, 'a')

pointsGroup = dbFile.create_group("points")
velocityGroup = dbFile.create_group("velocity")

pointsGroup.create_dataset("pointsY", data=pointsY)
pointsGroup.create_dataset("pointsZ", data=pointsY)

velocityGroup.create_dataset("uMean", data=uMean)
velocityGroup.create_dataset("uPrimeX", data=uPrimeX, compression="gzip")
velocityGroup.create_dataset("uPrimeY", data=uPrimeY, compression="gzip")
velocityGroup.create_dataset("uPrimeZ", data=uPrimeZ, compression="gzip")

pointsGroup.attrs["nPointsY"] = pointsY.shape[0]
pointsGroup.attrs["nPointsZ"] = pointsY.shape[1]
pointsGroup.attrs["nPoints"] = pointsY.size
