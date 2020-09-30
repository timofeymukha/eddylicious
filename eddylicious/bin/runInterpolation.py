# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

from __future__ import print_function
from __future__ import division
import os
import numpy as np
import argparse
import functools
from mpi4py import MPI
import h5py
from scipy.spatial import Delaunay
from eddylicious.helper_functions import config_to_dict, set_write_path
from eddylicious.readers.foamfile_readers import read_points_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_foamfile
from eddylicious.writers.ofnative_writers import write_points_to_ofnative, write_velocity_to_ofnative
from eddylicious.writers.hdf5_writers import write_points_to_hdf5, write_velocity_to_hdf5
from eddylicious.writers.vtk_writers import write_data_to_vtk
from eddylicious.generators.interpolation import interpolation_generate


def get_times(reader, readPath):
    """Read the time values associated with the precursor database."""

    # Grab the existing times and sort them
    if reader == "foamFile":
        times = os.listdir(readPath)
        times = np.sort(times)
    elif reader == "hdf5":
        # Set the readPath to the file itself
        readPath = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        times = readPath["velocity"]["times"][:]
        readPath.close()
    else:
        raise ValueError("Unknown reader: "+reader)

    return times


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

    if rank == 0:
        print("Producing database with "+str(size)+" time-steps.")


# SET UP GEOMETRY
# Read grid for the recycling plane

    times = get_times(reader, readPath)
    if reader == "foamFile":
        sampleSurfaceName = configDict["sampleSurfaceName"]
        sampleFunctionObjectName = configDict["sampleFunctionObjectName"]
        dataDir = os.path.join(readPath, "postProcessing", sampleFunctionObjectName)
        times = get_times(reader, dataDir)
        pointsReadPath = os.path.join(dataDir, times[0], sampleSurfaceName,
                                      "faceCentres")
        [pointsY, pointsZ] = read_points_foamfile(pointsReadPath)

    else:
        raise ValueError("Unsupported or unknown reader: "+reader)

    if rank == 0:
        print("Reading from database with "+str(len(times)) + " time-steps.")

    if size > len(times):
        raise ValueError("desired time-span is too large!")

    # Bounds for the points

    minYPrec = float(configDict["minYPrec"])
    maxYPrec = float(configDict["maxYPrec"])
    minZPrec = float(configDict["minZPrec"])
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
    minYInfl = float(configDict["minYInfl"])
    maxYInfl = float(configDict["maxYInfl"])
    minZInfl = float(configDict["minZInfl"])
    maxZInfl = float(configDict["maxZInfl"])

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
        dataDir = os.path.join(readPath, "postProcessing", sampleFunctionObjectName)
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
