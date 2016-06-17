from __future__ import print_function
from __future__ import division
import os
import numpy as np
import argparse
from mpi4py import MPI
import h5py
from eddylicious.generators.helper_functions import *
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_from_foamfile
from eddylicious.readers.hdf5_readers import read_points_from_hdf5
from eddylicious.readers.hdf5_readers import read_velocity_from_hdf5
from eddylicious.writers.ofnative_writers import write_points_to_ofnative
from eddylicious.writers.hdf5_writers import write_points_to_hdf5
from eddylicious.generators.lund_rescaling import lund_generate
from eddylicious.generators.lund_rescaling import lund_rescale_mean_velocity


def set_write_path(config):
    """Sets the writePath variable in concordance with the writer.

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


def get_umean_prec(reader, readPath, flip):
    """Reed the mean velocity profile of the precursor
       and the total number of points in the y direction.

    """
    if reader == "foamFile":
        uMeanTimes = os.listdir(os.path.join(readPath, "postProcessing",
                                             "collapsedFields"))
        uMean = np.genfromtxt(os.path.join(readPath, "postProcessing",
                                           "collapsedFields",
                                           uMeanTimes[-1],
                                           "UMean_X.xy"))
        uMeanX = uMean[:, 1]
        if uMean.shape[1] == 3:
            uMeanY = uMean[:, 2]
        else:
            uMeanY = np.zeros(uMeanX.shape)
    elif reader == "hdf5":
        readPath = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        uMeanX = readPath["velocity"]["uMeanX"][:]
        uMeanY = readPath["velocity"]["uMeanY"][:]
        readPath.close()
    else:
        raise ValueError("Unknown reader: "+reader)

    return uMeanX, uMeanY


def get_y_prec(reader, readPath):
    """Read the mean velocity profile of the precursor
       and the total number of points in the y direction.

    """
    if reader == "foamFile":
        uMeanTimes = os.listdir(os.path.join(readPath, "postProcessing",
                                             "collapsedFields"))
        y = np.genfromtxt(os.path.join(readPath, "postProcessing",
                                       "collapsedFields",
                                       uMeanTimes[-1],
                                       "UMean_X.xy"))[:, 0]
    elif reader == "hdf5":
        readPath = h5py.File(readPath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        y = readPath["points"]["pointsY"][:, 0]
        readPath.close()
    else:
        raise ValueError("Unknown reader: "+reader)

    return y


def compute_tbl_properties(y, uMean, nu, flip):
    """Compute various parameters of a TBL."""

    y = y[np.nonzero(y)]
    uMean = uMean[np.nonzero(uMean)]

    if not flip:
        theta = momentum_thickness(y, uMean)
        delta = delta_99(y, uMean)
        deltaStar = delta_star(y, uMean)
        uTau = np.sqrt(nu*uMean[0]/y[0])
        u0 = uMean[-1]
    else:
        theta = momentum_thickness(np.flipud(y), uMean[::-1])
        delta = delta_99(np.flipud(y), uMean[::-1])
        deltaStar = delta_star(np.flipud(y), uMean[::-1])
        uTau = np.sqrt(nu*uMean[-1]/y[-1])
        u0 = uMean[0]

    yPlus1 = np.min(y)*uTau/nu
    return theta, deltaStar, delta, uTau, u0, yPlus1


def compute_ninfl(etaInfl, etaPrec):
    """Compute the number of elements that constitute the TBL."""

    nInfl = 0
    for i in range(etaInfl.size):
        if etaInfl[0] < etaInfl[-1]:
            if etaInfl[i] <= np.max(etaPrec):
                nInfl += 1
        elif np.flipud(etaInfl)[i] <= np.max(etaPrec):
                nInfl += 1
    return nInfl


def print_tbl_properties(theta, deltaStar, delta, uTau, u0, nu,
                         yPlus1):
    """Print out various parameters of a TBL."""

    reTheta = theta*u0/nu
    reDelta99 = delta*u0/nu
    reDeltaStar = deltaStar*u0/nu
    reTau = uTau*delta/nu

    print("    Re_theta "+str(reTheta))
    print("    Re_delta* "+str(reDeltaStar))
    print("    Re_delta99 "+str(reDelta99))
    print("    Re_tau "+str(reTau))
    print(" ")
    print("    theta "+str(theta))
    print("    delta* "+str(deltaStar))
    print("    delta99 "+str(delta))
    print("    u_tau "+str(uTau))
    print("    U0 "+str(u0))
    print("    cf "+str(0.5*(uTau/u0)**2))
    print("    y+_1 "+str(yPlus1))


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


# Determine which half of the channel flow will be used
    if configDict["half"] == "top":
        flipPrec = True
    elif configDict["half"] == "bottom":
        flipPrec = False
    else:
        raise ValueError("half should be either 'top' of 'bottom'")

# Viscosity
    nuInfl = float(configDict["nuInflow"])
    nuPrec = float(configDict["nuPrecursor"])

# Freestream velocity
    u0Infl = float(configDict["Ue"])

    if "delta99" in configDict:
        deltaInfl = float(configDict["delta99"])
        cfInfl = 0.02*pow(nuInfl/(u0Infl*deltaInfl), 1.0/6)
    elif "theta" in configDict:
        thetaInfl = float(configDict["theta"])
        cfInfl = 0.013435*(thetaInfl*u0Infl/nuInfl - 373.83)**(-2/11)
    else:
        raise ValueError("The config file should provide delta99 or theta")

    # Friction velocity
    if configDict["uTauInflow"] == "compute":
        uTauInfl = u0Infl*np.sqrt(cfInfl/2)
    else:
        uTauInfl = float(configDict["uTauInflow"])


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

    # get the times in the precursor database
    times = get_times(reader, readPath)

    if rank == 0:
        print("Reading from database with "+str(len(times)) + " time-steps.")

    # Get the mean velocity for the precursor
    uMeanXPrec, uMeanYPrec = get_umean_prec(reader, readPath, flipPrec)

    yPrec = get_y_prec(reader, readPath)
    centerY = (yPrec[0] + yPrec[-1])/2
    centerY = 3
    totalPointsY = yPrec.size

    indY = np.argmin(abs(yPrec - centerY))

    if not flipPrec:
        uMeanXPrec = uMeanXPrec[:indY+1]
        uMeanYPrec = uMeanYPrec[:indY+1]
    else:
        uMeanXPrec = uMeanXPrec[indY+1:]
        uMeanYPrec = uMeanYPrec[indY+1:]

    nPointsY = uMeanXPrec.size

# SET UP GEOMETRY
# Read grid for the recycling plane

    times = get_times(reader, readPath)
    if reader == "foamFile":
        sampleSurfaceName = configDict["sampleSurfaceName"]
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
        pointsReadPath = os.path.join(dataDir, times[0], sampleSurfaceName,
                                      "faceCentres")
        if not flipPrec:
            [pointsY, pointsZ, yInd, zInd] = \
                read_points_from_foamfile(pointsReadPath, addValBot=yPrec[0],
                                          addValTop=yPrec[-1],
                                          excludeTop=totalPointsY-nPointsY,
                                          exchangeValTop=centerY) #!!
        else:
            [pointsY, pointsZ, yInd, zInd] = \
                read_points_from_foamfile(pointsReadPath, addValBot=yPrec[0],
                                          addValTop=yPrec[-1],
                                          excludeBot=totalPointsY-nPointsY,
                                          exchangeValBot=centerY)
    elif reader == "hdf5":

        if not flipPrec:
            [pointsY, pointsZ] = \
                read_points_from_hdf5(readPath,
                                      excludeTop=totalPointsY-nPointsY,
                                      exchangeValTop=centerY) #!!
        else:
            [pointsY, pointsZ] = \
                read_points_from_hdf5(readPath,
                                      excludeBot=totalPointsY-nPointsY,
                                      exchangeValBot=centerY)
    else:
        raise ValueError("Unknown reader: "+reader)

    [nPointsY, nPointsZ] = pointsY.shape

# Read grid for the inflow plane
    if inflowReader == "foamFile":
        [pointsYInfl, pointsZInfl] =\
            read_points_from_foamfile(os.path.join(inflowGeometryPath))[:2]
    else:
        raise ValueError("Unknown inflow reader: "+inflowReader)

    nPointsZInfl = pointsYInfl.shape[1]

# Get the grid points along y as 1d arrays for convenience
    yPrec = pointsY[:, 0]
    if not flipPrec:
        yPrec = np.abs(yPrec - yPrec[0])
    else:
        yPrec = np.abs(yPrec - yPrec[-1])

    yInfl = pointsYInfl[:, 0]
    yInfl = np.abs(yInfl - yOrigin)

# See whether the inflow TBL is up side down or not
    if yInfl[0] < yInfl[1]:
        flipInfl = False
    else:
        flipInfl = True

    # Outer scale coordinates
    [thetaPrec, deltaStarPrec, deltaPrec,
     uTauPrec, u0Prec, yPlus1Prec] = compute_tbl_properties(yPrec, uMeanXPrec,
                                                            nuPrec, flipPrec)
    if "delta99" in configDict:
        etaPrec = yPrec/deltaPrec
        etaInfl = yInfl/deltaInfl
        blending = blending_function(etaInfl)
        #blending = np.ones(blending.shape)
    else:
        etaPrec = yPrec/thetaPrec
        etaInfl = yInfl/thetaInfl
        blending = blending_function(etaInfl/8)

# Inner scale coordinates
    gamma = uTauInfl/uTauPrec
    yPlusPrec = yPrec*uTauPrec/nuPrec
    yPlusInfl = yInfl*uTauInfl/nuInfl

# Points containing the boundary layer at the inflow plane
    nInfl = compute_ninfl(etaInfl, etaPrec)

# Create the reader functions
    if reader == "foamFile":
        dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
        if not flipPrec:
            readerFunc = read_velocity_from_foamfile(
                            dataDir,
                            sampleSurfaceName,
                            nPointsZ, yInd, zInd,
                            addValBot=(0, 0, 0), addValTop=(0, 0, 0),
                            excludeTop=totalPointsY-nPointsY,
                            interpValTop=True)
        else:
            readerFunc = read_velocity_from_foamfile(
                            dataDir,
                            sampleSurfaceName,
                            nPointsZ, yInd, zInd,
                            addValBot=(0, 0, 0), addValTop=(0, 0, 0),
                            excludeBot=totalPointsY-nPointsY,
                            interpValBot=True)
    elif reader == "hdf5":
        if not flipPrec:
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


    uMeanXInfl, uMeanYInfl = lund_rescale_mean_velocity(etaPrec, yPlusPrec,
                                                        uMeanXPrec, uMeanYPrec,
                                                        nInfl,
                                                        etaInfl, yPlusInfl,
                                                        nPointsZInfl,
                                                        u0Infl, u0Prec,
                                                        gamma,
                                                        blending)

    if rank == 0:
        print("Precursor properties:")
        print_tbl_properties(thetaPrec, deltaStarPrec, deltaPrec, uTauPrec,
                             u0Prec, nuPrec, yPlus1Prec)

# Generate the inflow fields
    if rank == 0:
        print("Generating the inflow fields.")

    comm.Barrier()

    lund_generate(readerFunc,
                  writer, writePath,
                  dt, t0, tEnd, timePrecision,
                  uMeanXPrec, uMeanXInfl,
                  uMeanYPrec, uMeanYInfl,
                  etaPrec, yPlusPrec, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  nInfl, gamma,
                  times, blending)

    if rank == 0:
        print("Process 0 done, waiting for the others...")

    comm.Barrier()

    if rank == 0:
        print("Done\n")

    [thetaInfl, deltaStarInfl, deltaInfl,
     uTauInfl, u0Infl, yPlus1Infl] = compute_tbl_properties(yInfl,
                                                            uMeanXInfl[:, 0],
                                                            nuInfl,
                                                            flipInfl)

    if rank == 0:
        print("Inflow boundary properties")
        print_tbl_properties(thetaInfl, deltaStarInfl, deltaInfl, uTauInfl,
                             u0Infl, nuInfl, yPlus1Infl)

if __name__ == "__main__":
    main()
