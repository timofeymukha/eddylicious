#!/pdc/vol/anaconda/2.1/py27/bin/python
from __future__ import print_function
from __future__ import division
import os
import numpy as np
import argparse
from mpi4py import MPI
import h5py as h5py
from eddylicious.generators.helper_functions import delta_99
from eddylicious.generators.helper_functions import theta
from eddylicious.generators.helper_functions import delta_star
from eddylicious.generators.helper_functions import blending_function
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_from_foamfile
from eddylicious.readers.hdf5_readers import read_points_from_hdf5
from eddylicious.readers.hdf5_readers import read_velocity_from_hdf5
from eddylicious.writers.tvmfv_writers import write_points_to_tvmfv
from eddylicious.writers.hdf5_writers import write_points_to_hdf5
from eddylicious.generators.lund_rescaling import lund_generate
from eddylicious.generators.lund_rescaling import lund_rescale_mean_velocity


def config_to_dict(configFile):
    # Read the config file into a dictionary
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
    inflowReadPath = configDict["inflowGeometryPath"]
    writePath = configDict["writePath"]

    sampleSurfaceName = configDict["sampleSurfaceName"]
    inletPatchName = configDict["inletPatchName"]

    reader = configDict["reader"]
    inflowReader = configDict["inflowGeometryReader"]
    writer = configDict["writer"]

    if writer == "hdf5":
        hdf5FileName = configDict["hdf5FileName"]

# Determine which half of the channel flow will be used
    if configDict["half"] == "top":
        flip = True
    elif configDict["half"] == "bottom":
        flip = False
    else:
        raise ValueError("half should be either 'top' of 'bottom'")

# Viscosity
    nuInfl = float(configDict["nuInflow"])
    nuPrec = float(configDict["nuPrecursor"])

# The outer scale
    if "delta99" in configDict:
        deltaInfl = float(configDict["delta99"])
        blendingFunction = blending_function
    elif "theta" in configDict:
        thetaInfl = float(configDict["theta"])
    else:
        raise ValueError("The config file should provide delta99 or theta")

# Freestream velocity
    u0Infl = float(configDict["Ue"])

# Friction velocities
    if "delta99" in configDict:
        reDelta99Infl = u0Infl*deltaInfl/nuInfl
        cfInfl = 0.02*pow(1.0/reDelta99Infl, 1.0/6)
    elif "theta" in configDict:
        reThetaInfl = u0Infl*thetaInfl/nuInfl

    if configDict["uTauInflow"] == "compute":
# TODO: TREAT case when theta is provided
        uTauInfl = u0Infl*np.sqrt(cfInfl/2)
    else:
        uTauInfl = float(configDict["uTauInflow"])

    uTauPrec = float(configDict["uTauPrecursor"])

    gamma = uTauInfl/uTauPrec

# Shift in coordinates
    xOrigin = float(configDict["xOrigin"])
    yOrigin = float(configDict["yOrigin"])

# Time-step and initial time for the writer
    dt = float(configDict["dt"])
    t0 = float(configDict["t0"])
    tEnd = float(configDict["tEnd"])
    timePrecision = int(configDict["tPrecision"])
    size = int((tEnd-t0)/dt+1)

    if rank == 0:
        print("Producing database with "+str(size)+" time-steps.")

# SET UP THE READER

    if reader == "foamFile":
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
    elif reader == "hdf5":
        dataDir = readPath
    else:
        raise ValueError("Unknown reader: "+reader)

# Grab the existing times and sort them
    if reader == "foamFile":
        times = os.listdir(dataDir)
        times = np.sort(times)
    elif reader == "hdf5":
        dbFile = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        times = dbFile["velocity"]["times"]
    else:
        raise ValueError("Unknown reader: "+reader)

    if rank == 0:
        print("Reading from database with "+str(len(times)) + " time-steps.")

# SET UP MEAN PROFILE
    if reader == "foamFile":
        uMeanTimes = os.listdir(os.path.join(readPath, "postProcessing",
                                             "collapsedFields"))
        uMean = np.genfromtxt(os.path.join(readPath, "postProcessing",
                                           "collapsedFields",
                                           uMeanTimes[-1],
                                           "UMean_X.xy"))[:, 1]
        uMean = np.append(np.zeros((1, 1)), uMean)
        uMean = np.append(uMean, np.zeros((1, 1)))
    elif reader == "hdf5":
        uMean = dbFile["velocity"]["uMean"]
    else:
        raise ValueError("Unknown reader: "+reader)

    totalPointsY = uMean.size

    if not flip:
        uMean = uMean[:int(totalPointsY*0.5)]
    else:
        uMean = uMean[int(totalPointsY*0.5)+1:]

    u0Prec = np.max(uMean)

    nPointsY = uMean.size

# SET UP GEOMETRY
# Read grid for the recycling plane
    if reader == "foamFile":
        pointsReadPath = os.path.join(dataDir, times[0], sampleSurfaceName,
                                      "faceCentres")
        if not flip:
            [pointsY, pointsZ, yInd, zInd] = \
                read_points_from_foamfile(pointsReadPath, addValBot=0,
                                          addValTop=2,
                                          excludeTop=totalPointsY-nPointsY,
                                          exchangeValTop=1.0)
        else:
            [pointsY, pointsZ, yInd, zInd] = \
                read_points_from_foamfile(pointsReadPath, addValBot=0,
                                          addValTop=2,
                                          excludeBot=totalPointsY-nPointsY,
                                          exchangeValBot=1.0)
    elif reader == "hdf5":

        if not flip:
            [pointsY, pointsZ] = \
                read_points_from_hdf5(readPath,
                                      excludeTop=totalPointsY-nPointsY,
                                      exchangeValTop=1.0)
        else:
            [pointsY, pointsZ] = \
                read_points_from_hdf5(readPath,
                                      excludeBot=totalPointsY-nPointsY,
                                      exchangeValBot=1.0)
    else:
        raise ValueError("Unknown reader: "+reader)

    [nPointsY, nPointsZ] = pointsY.shape

# Read grid for the inflow plane
    if inflowReader == "foamFile":
        [pointsYInfl, pointsZInfl, yIndInfl, zIndInfl] =\
            read_points_from_foamfile(os.path.join(inflowReadPath))
    else:
        raise ValueError("Unknown inflow reader: "+inflowReader)

    nPointsZInfl = pointsYInfl.shape[1]

# Check if Lz is the same
    if (not pointsZInfl[0, -1] == pointsZ[0, -1]) and (rank == 0):
        print("Warning: the lengths of the domains in the spanwise direction"
              "is not equal")


# Get the grid points along y as 1d arrays for convenience
    yPrec = pointsY[:, 0]
    yInfl = pointsYInfl[:, 0]

# Outer scale coordinates
    if not flip:
        yOriginPrec = 0
        deltaPrec = delta_99(yPrec, uMean)
        thetaPrec = theta(yPrec, uMean)
    else:
    # TODO should be generic with respect to precursor channel height
        yOriginPrec = 2
        deltaPrec = delta_99(np.flipud(np.abs(yPrec - yOrigin)),
                             np.flipud(uMean))
        thetaPrec = theta(np.flipud(np.abs(yPrec - yOrigin)),
                          np.flipud(uMean))

    if "delta99" in configDict:
        etaPrec = np.abs(yPrec - yOriginPrec)/deltaPrec
        etaInfl = np.abs(yInfl - yOrigin)/deltaInfl
    elif "theta" in configDict:
        etaPrec = np.abs(yPrec - yOriginPrec)/thetaPrec
        etaInfl = np.abs(yInfl - yOrigin)/thetaInfl


# Inner scale coordinates
    yPlusPrec = np.abs(yPrec - yOriginPrec)*uTauPrec/nuPrec
    yPlusInfl = np.abs(yInfl - yOrigin)*uTauInfl/nuInfl

# Points containing the boundary layer at the inflow plane
    nInfl = 0
    for i in range(etaInfl.size):
        if etaInfl[0] < etaInfl[-1]:
            if etaInfl[i] <= np.max(etaPrec):
                nInfl += 1
        elif np.flipud(etaInfl)[i] <= np.max(etaPrec):
                nInfl += 1

# Points where inner scaling will be used
    nInner = 0
    for i in range(etaInfl.size):
        if etaInfl[i] <= 0.7:
            nInner += 1

# Create the reader functions
    if reader == "foamFile":
        if not flip:
            readerFunc = read_velocity_from_foamfile(
                            dataDir,
                            sampleSurfaceName,
                            nPointsZ, yInd, zInd,
                            addValBot=0, addValTop=0,
                            excludeTop=totalPointsY-nPointsY,
                            interpValTop=True)
        else:
            readerFunc = read_velocity_from_foamfile(
                            dataDir,
                            sampleSurfaceName,
                            nPointsZ, yInd, zInd,
                            addValBot=0, addValTop=0,
                            excludeBot=totalPointsY-nPointsY,
                            interpValBot=True)
    elif reader == "hdf5":
        if not flip:
            readerFunc = read_velocity_from_hdf5(
                            readPath,
                            excludeTop=totalPointsY-nPointsY,
                            interpValTop=True)
        else:
            readerFunc = read_velocity_from_hdf5(
                            readPath,
                            excludeBot=totalPointsY-nPointsY,
                            interpValBot=True)
    else:
        raise ValueError("Unknown reader: "+reader)


# SET UP WRITERS

# Write points and modify writePath appropriately
    if writer == "tvmfv":
        writePath = os.path.join(writePath, "constant", "boundaryData",
                                 inletPatchName)
        if rank == 0:
            if not os.path.exists(writePath):
                os.makedirs(writePath)

        write_points_to_tvmfv(os.path.join(writePath, "points"), pointsYInfl,
                              pointsZInfl, xOrigin)

    elif writer == "hdf5":
        writePath = os.path.join(writePath, hdf5FileName)
# If the hdf5 file exists, delete it.
        if rank == 0 and os.path.isfile(writePath):
            print("HDF5 database already exists. It it will be overwritten.")
            os.remove(writePath)

        if rank == 0:
            write_points_to_hdf5(writePath, pointsYInfl, pointsZInfl, xOrigin)

# We change the writePath to be the hdf5 file itself
        writePath = h5py.File(writePath, 'a', driver='mpio',
                              comm=MPI.COMM_WORLD)
        writePath.create_dataset("time", data=t0*np.ones((size, 1)))
        writePath.create_dataset("velocity", (size, pointsZInfl.size, 3),
                                 dtype=np.float64)
    else:
        raise ValueError("Unknown writer: "+writer)

    uMeanInfl = lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMean,
                                           nInfl, nInner,
                                           etaInfl, yPlusInfl, nPointsZInfl,
                                           u0Infl, u0Prec, gamma,
                                           blendingFunction)

    if not flip:
        if "delta99" in configDict:
            thetaInfl = theta(yInfl, uMeanInfl[:, 0])
        else:
            deltaInfl = delta_99(yInfl, uMeanInfl[:, 0])
        deltaStarInfl = delta_star(yInfl, uMeanInfl[:, 0])
    else:
        if "delta99" in configDict:
            thetaInfl = theta(np.flipud(np.abs(yInfl - yOrigin)),
                              uMeanInfl[::-1, 0])
        else:
            deltaInfl = delta_99(np.flipud(np.abs(yInfl - yOrigin)),
                                 uMeanInfl[::-1, 0])
        deltaStarInfl = delta_star(np.flipud(np.abs(yInfl - yOrigin)),
                                   uMeanInfl[::-1, 0])

    if "delta99" in configDict:
        reThetaInfl = thetaInfl*u0Infl/nuInfl
    elif "theta" in configDict:
        reDelta99Infl = deltaInfl*u0Infl/nuInfl

    reDeltaStarInfl = deltaStarInfl*u0Infl/nuInfl
    reTauInfl = uTauInfl*deltaInfl/nuInfl

# Generate the inflow fields
    if rank == 0:
        print("Generating the inflow fields.")

    lund_generate(readerFunc,
                  writer, writePath,
                  dt, t0, tEnd, timePrecision,
                  uMean, uMeanInfl,
                  etaPrec, yPlusPrec, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  nInfl, nInner, gamma,
                  times, blendingFunction)

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()
    if rank == 0:
        print("Done\n")

    if rank == 0:
        print("Inflow boundary properties:")
        print("    Re_theta "+str(reThetaInfl))
        print("    Re_delta* "+str(reDeltaStarInfl))
        print("    Re_delta99 "+str(reDelta99Infl))
        print("    Re_tau "+str(reTauInfl))
        print(" ")
        print("    theta "+str(thetaInfl))
        print("    delta* "+str(deltaStarInfl))
        print("    delta99 "+str(deltaInfl))

if __name__ == "__main__":
    main()
