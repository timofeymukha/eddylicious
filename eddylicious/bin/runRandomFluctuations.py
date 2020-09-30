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
from eddylicious.helper_functions import config_to_dict, set_write_path
from eddylicious.readers.foamfile_readers import read_points_foamfile
from eddylicious.writers.ofnative_writers import write_points_to_ofnative, write_velocity_to_ofnative
from eddylicious.writers.hdf5_writers import write_points_to_hdf5, write_velocity_to_hdf5
from eddylicious.generators.random_fluctuations import random_fluctuations_generate


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

# Parse the command-line arguments
    parser = argparse.ArgumentParser(
            description="A script for generating inflow \
                         random velocity fields with predetermined \
                         first- and second-order statistics.")

    parser.add_argument('--config',
                        type=str,
                        help='The config file',
                        required=True)

    args = parser.parse_args()

# PARSE THE CONFIG
    configFile = open(args.config, mode='r')
    configDict = config_to_dict(configFile)


# Readers and writers
    uMeanPath = configDict["uMeanPath"]
    reynoldsStressPath = configDict["reynoldsStressPath"]
    inflowGeometryPath = configDict["inflowGeometryPath"]
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

# Read grid for the inflow plane
    if inflowReader == "foamFile":
        [pointsYInfl, pointsZInfl] = \
            read_points_foamfile(os.path.join(inflowGeometryPath))
    else:
        raise ValueError("Unknown inflow reader: "+inflowReader)

    if rank == 0:
        print("\nTotal number of points of the inflow patch", pointsYInfl.shape[0])

    # Get the write path appropriate for the reader
    writePath = set_write_path(configDict)

    # Create the writer functions
    if writer == "ofnative":
        if rank == 0:
            write_points_to_ofnative(os.path.join(writePath, "points"),
                                     pointsYInfl, pointsZInfl, xOrigin)
        writerFunc = functools.partial(write_velocity_to_ofnative, writePath)
    elif writer == "hdf5":
        writePath.create_dataset("time", data=t0*np.ones((size, 1)))
        writePath.create_dataset("velocity", (size, pointsZInfl.size, 3),
                                 dtype=np.float64)
        write_points_to_hdf5(writePath, pointsYInfl, pointsZInfl, xOrigin)
        writerFunc = functools.partial(write_velocity_to_hdf5, writePath)
    else:
        raise ValueError("Unknown writer: "+writer)

# Read the mean velocity and reynolds stress distributions
    uMean = np.genfromtxt(uMeanPath)
    reynoldsStress = np.genfromtxt(reynoldsStressPath)

# Generate the inflow fields
    if rank == 0:
        print("Generating")

    comm.Barrier()

    random_fluctuations_generate(writerFunc,
                                 np.column_stack((pointsYInfl, pointsZInfl)),
                                 uMean,
                                 reynoldsStress,
                                 dt, t0, tEnd, timePrecision)

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()

    if rank == 0:
        print("Done\n")


if __name__ == "__main__":
    main()
