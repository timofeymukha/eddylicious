# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

from __future__ import print_function
from __future__ import division
import os
import numpy as np
import argparse
from mpi4py import MPI
import h5py
from scipy.spatial import Delaunay
from eddylicious.generators.helper_functions import *
from eddylicious.readers.foamfile_readers import read_points_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_foamfile
from eddylicious.readers.hdf5_readers import read_structured_points_hdf5
from eddylicious.readers.hdf5_readers import read_structured_velocity_hdf5
from eddylicious.writers.ofnative_writers import write_points_to_ofnative
from eddylicious.writers.hdf5_writers import write_points_to_hdf5
from eddylicious.generators.interpolation import interpolation_generate


def set_write_path(config):
    """Set the writePath variable in concordance with the writer.

    For the ofnative writer: the path to constant/boundaryData directory.
    For the hdf5 writer: the hdf5 file itself.

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    writer = config["writer"]
    writePath = config["writePath"]

    if writer == "ofnative":
        inletPatchName = config["inletPatchName"]
        writePath = os.path.join(writePath, "constant", "boundaryData",
                                 inletPatchName)
        if rank == 0:
            if not os.path.exists(writePath):
                os.makedirs(writePath)

    elif writer == "hdf5":
        writePath = os.path.join(writePath, config["hdf5FileName"])
        # If the hdf5 file exists, delete it.
        if rank == 0 and os.path.isfile(writePath):
            print("HDF5 database already exists. It it will be overwritten.")
            os.remove(writePath)

        # We change the writePath to be the hdf5 file itself
        writePath = h5py.File(writePath, 'a', driver='mpio',
                              comm=MPI.COMM_WORLD)
    else:
        raise ValueError("Unknown writer: "+writer)

    return writePath


def get_times(reader, readPath):
    """Read the time values associated with the precursor database."""

    # Grab the existing times and sort them
    if reader == "foamFile":
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
        times = os.listdir(dataDir)
        times = np.sort(times)
    elif reader == "hdf5":
        # Set the readPath to the file itself
        readPath = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        times = readPath["velocity"]["times"][:]
        readPath.close()
    else:
        raise ValueError("Unknown reader: "+reader)

    return times


def config_to_dict(configFile):
    """Parse a config file to a dictionary."""

    configDict = {}

    for line in configFile:
        if (line[0] == '#') or (line == '\n'):
            continue
        configDict[line.split()[0]] = line.split()[1]

    return configDict


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

# Parse the command-line arguments
    parser = argparse.ArgumentParser(
            description="A script for generating inflow \
                            velocity fields using Lund et al's rescaling.")

    parser.add_argument('--config',
                        type=str,
                        help='The config file',
                        required=True)

    args = parser.parse_args()

# PARSE THE CONFIG
    configFile = open(args.config, mode='r')
    configDict = config_to_dict(configFile)


# Readers and writers
    readPath = configDict["readPath"]
    inflowGeometryPath = configDict["inflowGeometryPath"]

    reader = configDict["reader"]
    inflowReader = configDict["inflowGeometryReader"]
    writer = configDict["writer"]



# Shift x in coordinates
    if "xOrigin" in configDict:
        xOrigin = float(configDict["xOrigin"])
    else:
        xOrigin = 0.0

# Time-step and initial time for the writer
    dt = float(configDict["dt"])
    t0 = float(configDict["t0"])
    tEnd = float(configDict["tEnd"])
    timePrecision = int(configDict["tPrecision"])
    size = int((tEnd-t0)/dt+1)

    # get the times in the precursor database
    times = get_times(reader, readPath)

    if rank == 0:
        print("Producing database with "+str(size)+" time-steps.")
        print("Reading from database with "+str(len(times)) + " time-steps.")

    if size > len(times):
        raise ValueError("desired time-span is too large!")

# SET UP GEOMETRY
# Read grid for the recycling plane

    times = get_times(reader, readPath)
    if reader == "foamFile":
        sampleSurfaceName = configDict["sampleSurfaceName"]
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
        pointsReadPath = os.path.join(dataDir, times[0], sampleSurfaceName,
                                      "faceCentres")
        [pointsY, pointsZ] = read_points_foamfile(pointsReadPath)

    else:
        raise ValueError("Unsupported or unknown reader: "+reader)

    # Bounds for the points

    if "minYPrec" in configDict:
        minYPrec = float(configDict["minYPrec"])
    else:
        minYPrec = pointsY.min()

    if "maxYPrec" in configDict:
        maxYPrec = float(configDict["maxYPrec"])
    else:
        maxYPrec = pointsY.max()

    if "minZPrec" in configDict:
        minZPrec = float(configDict["minZPrec"])
    else:
        minZPrec = pointsZ.min()


    if "maxZPrec" in configDict:
        maxZPrec = float(configDict["maxZPrec"])
    else:
        maxZPrec = pointsZ.max()

    if rank == 0:
        print("\nTotal number of source points,", pointsY.shape[0]) 
        print("Filtering source points with bounding box, y:[", minYPrec,
              maxYPrec, "] z:[", minZPrec, maxZPrec, "]")

    idxPrec = np.where((pointsY >= minYPrec) & (pointsY <= maxYPrec) &
                       (pointsZ >= minZPrec) & (pointsY <= maxZPrec))

    # Normalize to unit square
    pointsY = (pointsY[idxPrec] - minYPrec)/(maxYPrec - minYPrec)
    pointsZ = (pointsZ[idxPrec] - minZPrec)/(maxZPrec - minZPrec)

    if rank == 0:
        print("Source points after filtering:", pointsY.shape[0])

    triangulation = Delaunay(np.column_stack((pointsY, pointsZ)))

# Read grid for the inflow plane
    if inflowReader == "foamFile":
        [pointsYInfl, pointsZInfl] = \
            read_points_foamfile(os.path.join(inflowGeometryPath))
    else:
        raise ValueError("Unknown inflow reader: "+inflowReader)

    # Bound the points
    if "minYInfl" in configDict:
        minYInfl = float(configDict["minYInfl"])
    else:
        minYInfl = pointsYInfl.min()

    if "maxYInfl" in configDict:
        maxYInfl = float(configDict["maxYInfl"])
    else:
        maxYInfl = pointsYInfl.max()

    if "minZInfl" in configDict:
        minZInfl = float(configDict["minZInfl"])
    else:
        minZInfl = pointsZInfl.min()

    if "maxZInfl" in configDict:
        maxZInfl = float(configDict["maxZInfl"])
    else:
        maxZInfl = pointsZInfl.max()

    if rank == 0:
        print("\nTotal number of target points,", pointsYInfl.shape[0]) 
        print("Filtering target points with bounding box, y:[", minYInfl,
              maxYInfl, "] z:[", minZInfl, maxZInfl, "]")

    # Filter the points to those inside rectangle
    idxInfl = np.where((pointsYInfl >= minYInfl) & (pointsYInfl <= maxYInfl) &
                       (pointsZInfl >= minZInfl) & (pointsZInfl <= maxZInfl))


    pointsYInfl = pointsYInfl[idxInfl]
    pointsZInfl = pointsZInfl[idxInfl]

    if rank == 0:
        print("Target points after filtering:", pointsYInfl.shape[0])

# Create the reader functions
    if reader == "foamFile":
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
        readerFunc = read_velocity_foamfile(dataDir, sampleSurfaceName)
    else:
        raise ValueError("Unsupported or unknown reader: "+reader)

    # Get the write path appropriate for the reader
    writePath = set_write_path(configDict)

    if writer == "ofnative":
        if rank == 0:
            write_points_to_ofnative(os.path.join(writePath, "points"),
                                     pointsYInfl, pointsZInfl, xOrigin)
    elif writer == "hdf5":
        writePath.create_dataset("time", data=t0*np.ones((size, 1)))
        writePath.create_dataset("velocity", (size, pointsZInfl.size, 3),
                                 dtype=np.float64)
        write_points_to_hdf5(writePath, pointsYInfl, pointsZInfl, xOrigin)

    # Transform inflow points to square
    pointsYInfl = (pointsYInfl - minYInfl)/(maxYInfl - minYInfl)
    pointsZInfl = (pointsZInfl - minZInfl)/(maxZInfl - minZInfl)


# Generate the inflow fields
    if rank == 0:
        print("Interpolating")

    comm.Barrier()

    interpolation_generate(readerFunc,
                           writer, writePath,
                           dt, t0, tEnd, timePrecision,
                           triangulation,
                           np.column_stack((pointsYInfl,pointsZInfl)),
                           idxPrec,
                           times)

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()

    if rank == 0:
        print("Done\n")

if __name__ == "__main__":
    main()
