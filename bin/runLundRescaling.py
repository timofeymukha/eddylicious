#!/pdc/vol/anaconda/2.1/py27/bin/python
import os
import numpy as np
from sys import exit
import argparse
from mpi4py import MPI
import h5py as h5py
from eddylicious.generators.helper_functions import delta_99
from eddylicious.generators.helper_functions import theta
from eddylicious.generators.helper_functions import delta_star
from eddylicious.readers.foamfile_readers import read_points_from_foamfile
from eddylicious.readers.foamfile_readers import read_velocity_from_foamfile
from eddylicious.readers.hdf5_readers import read_points_from_hdf5
from eddylicious.readers.hdf5_readers import read_velocity_from_hdf5
from eddylicious.writers.tvmfv_writers import write_points_to_tvmfv
from eddylicious.writers.hdf5_writers import write_points_to_hdf5
from eddylicious.generators.lund_rescaling import lund_generate
from eddylicious.generators.lund_rescaling import lund_rescale_mean_velocity

# Get processor rank
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

# Read the config file into a dictionary
configDict = {}

configFile = open(args.config, mode='r')

for line in configFile:
    if (line[0] == '#') or (line == '\n'):
        continue
    configDict[line.split()[0]] = line.split()[1]

# Read the readers and writers
readPath = configDict["readPath"]
inflowReadPath = configDict["inflowReadPath"]
writePath = configDict["writePath"]

sampleSurfaceName = configDict["sampleSurfaceName"]
inletPatchName = configDict["inletPatchName"]

reader = configDict["reader"]
inflowReader = configDict["inflowReader"]
writer = configDict["writer"]

hdf5FileName = configDict["hdf5FileName"]

# Determine which half of the channel flow will be used
if configDict["half"] == "top":
    flip = True
elif configDict["half"] == "bottom":
    flip = False
else:
    raise ValueError("half should be either 'top' of 'bottom'")

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
    print "Reading from database with ", len(times), " time-steps."

# Get the mean profile
if reader == "foamFile":
    uMeanTimes = os.listdir(os.path.join(readPath, "postProcessing",
                                         "collapsedFields"))
    uMean =  np.genfromtxt(os.path.join(readPath, "postProcessing",
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

if flip == False:
    uMean = uMean[:int(totalPointsY*0.5)]
else:
    uMean = uMean[int(totalPointsY*0.5)+1:]

nPointsY = uMean.size

# Read grid for the recycling plane
if reader == "foamFile":
    pointsReadPath = os.path.join(dataDir, times[0], sampleSurfaceName,
                                  "faceCentres")

    if not flip:
        [pointsY, pointsZ, yInd, zInd] = \
            read_points_from_foamfile(pointsReadPath, addValBot=0, addValTop=2,
                                    excludeTop=totalPointsY-nPointsY,
                                    exchangeValTop=1.0)
    else:
        [pointsY, pointsZ, yInd, zInd] = \
            read_points_from_foamfile(pointsReadPath, addValBot=0, addValTop=2,
                                    excludeBot=totalPointsY-nPointsY,
                                    exchangeValBot=1.0)

elif reader == "hdf5":

    if not flip:
        [pointsY, pointsZ] = \
            read_points_from_hdf5(readPath,  addValBot=0, addValTop=2,
                                  excludeTop=totalPointsY-nPointsY,
                                  exchangeValTop=1.0)
    else:
        [pointsY, pointsZ, yInd, zInd] = \
            read_points_from_foamfile(readPath, addValBot=0, addValTop=2,
                                      excludeBot=totalPointsY-nPointsY,
                                      exchangeValBot=1.0)
else:
    raise ValueError("Unknown reader: "+reader)


[nPointsY, nPointsZ] = pointsY.shape

# Read grid for the inflow plane
if inflowReader == "foamFile":
    [pointsYInfl, pointsZInfl, yIndInfl, zIndInfl] = read_points_from_foamfile(
        os.path.join(inflowReadPath, inletPatchName, "faceCentres"))

else:
    raise ValueError("Unknown inflow reader: "+inflowReader)

[nPointsYInfl, nPointsZInfl] = pointsYInfl.shape

# Define and compute the parameters of the boundary layer

# Viscosity
nuInfl = float(configDict["nuInflow"])
nuPrec = float(configDict["nuPrecursor"])


# Boundary layer thickneses
deltaInfl = float(configDict["delta99"])

# Freestream velocity, and centerline velocity
U0 = np.max(uMean)
Ue = float(configDict["Ue"])

# Reynolds number based on delta
ReDelta = Ue*deltaInfl/nuInfl

# Cf at the inflow, to get u_tau
cfInfl = 0.02*pow(1.0/ReDelta, 1.0/6)

# Friction velocities
if configDict["uTauInflow"] == "compute":
    uTauInfl = Ue*np.sqrt(cfInfl/2)
else:
    uTauInfl = float(configDict["uTauInflow"])

uTauPrec = float(configDict["uTauPrecursor"])

ReTauInfl = uTauInfl*deltaInfl/nuInfl

# gamma, the ratio of friction velocities
gamma = uTauInfl/uTauPrec

# Shift in coordinates
xOrigin = float(configDict["xOrigin"])
yOrigin = float(configDict["yOrigin"])

# Time-step and initial time for the writer
dt = float(configDict["dt"])
t0 = float(configDict["t0"])
tEnd = float(configDict["tEnd"])
timePrecision = int(configDict["tPrecision"])

# Get the grid points along y as 1d arrays for convenience
yPrec = pointsY[:, 0]
yInfl = pointsYInfl[:, 0]

# Outer scale coordinates
if flip == False:
    yOriginPrec = 0
    deltaPrec = delta_99(yPrec, uMean)
else:
    yOriginPrec = 2
    deltaPrec = delta_99(np.flipud(np.abs(yPrec-2)), np.flipud(uMean))

ReTauPrec = uTauPrec*deltaPrec/nuPrec

etaPrec = np.abs(yPrec - yOriginPrec)/deltaPrec
etaInfl = np.abs(yInfl- yOrigin)/deltaInfl

# Inner scale coordinates
yPlusPrec = np.abs(yPrec - yOriginPrec)*uTauPrec/nuPrec
yPlusInfl = np.abs(yInfl - yOrigin)*uTauInfl/nuInfl

# Points containing the boundary layer at the inflow plane
nInfl = 0
for i in xrange(etaInfl.size):  
    if etaInfl[0] < etaInfl[-1]:
        if etaInfl[i] <= np.max(etaPrec):
            nInfl += 1
    elif np.flipud(etaInfl)[i] <= np.max(etaPrec):
            nInfl += 1

# Points where inner scaling will be used
nInner = 0
for i in xrange(etaInfl.size):
    if etaInfl[i] <= 0.7:
        nInner += 1

if ReTauInfl > ReTauPrec and rank == 0:
    print "WARNING: Re_tau in the precursor is lower than in the desired TBL"

size = int((tEnd-t0)/dt+1)

if rank == 0:
    print "Producing database with", size, "time-steps."

# Write points and modify writePath appropriately
if writer == "tvmfv":
    if rank == 0:
        write_points_to_tvmfv(os.path.join(writePath, "constant",
                                           "boundaryData", inletPatchName,
                                           "points"),
                              pointsYInfl, pointsZInfl, xOrigin)
    writePath = os.path.join(writePath, "constant", "boundaryData",
                             inletPatchName)
elif writer == "hdf5":
    writePath = os.path.join(writePath, hdf5FileName)
    # If the hdf5 file exists, delete it.
    if rank == 0:
        if os.path.isfile(writePath):
            print "HDF5 database already exsists. It it will be overwritten."
            os.remove(writePath)
        write_points_to_hdf5(writePath, pointsYInfl, pointsZInfl, xOrigin)

    # We change the writePath to be the hdf5 file itself
    writePath = h5py.File(writePath, 'a', driver='mpio', comm=MPI.COMM_WORLD)
    writePath.create_dataset("time", data=t0*np.ones((size, 1)))
    writePath.create_dataset("velocity", (size, pointsZInfl.size, 3),
                             dtype=np.float)
else:
    raise ValueError("Unknown writer: "+writer)

# Create the reader functions
if reader == "foamFile":
    if not flip:
        readerFunc = read_velocity_from_foamfile(dataDir, sampleSurfaceName,
                                                 nPointsZ, yInd, zInd,
                                                 addValBot=0, addValTop=0,
                                                 excludeTop=totalPointsY-nPointsY,
                                                 interpValTop=True)
    else:
        readerFunc = read_velocity_from_foamfile(dataDir, sampleSurfaceName,
                                                 nPointsZ, yInd, zInd,
                                                 addValBot=0, addValTop=0,
                                                 excludeBot=totalPointsY-nPointsY,
                                                 interpValBot=True)
elif reader == "hdf5":
    if not flip:
        readerFunc = read_velocity_from_hdf5(readPath,
                                             addValBot=0, addValTop=0,
                                             excludeTop=totalPointsY-nPointsY,
                                             interpValTop=True)
    else:
        readerFunc = read_velocity_from_hdf5(readPath,
                                             addValBot=0, addValTop=0,
                                             excludeBot=totalPointsY-nPointsY,
                                             interpValBot=True)
else:
    raise ValueError("Unknown reader: "+reader)

uMeanInfl = lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMean,
                                       nInfl, nInner,
                                       etaInfl, yPlusInfl, nPointsZInfl,
                                       Ue, U0, gamma)

ReThetaInfl = theta(yInfl, uMeanInfl[:, 0])*Ue/nuInfl
ReDeltaStarInfl = delta_star(yInfl, uMeanInfl[:, 0])*Ue/nuInfl

# Generate the inflow fields
if rank == 0:
    print "Generating the inflow fields."

lund_generate(readerFunc,
              writer, writePath,
              dt, t0, tEnd, timePrecision,
              uMean, uMeanInfl,
              etaPrec, yPlusPrec, pointsZ,
              etaInfl, yPlusInfl, pointsZInfl,
              nInfl, nInner, gamma,
              times)
if rank == 0:
    print "Process 0 done, waiting for the others..."
