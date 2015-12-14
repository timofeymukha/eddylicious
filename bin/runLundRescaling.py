import os
import numpy as np
from sys import exit
import argparse
from eddylicous.generators.helper_functions import delta_99
from eddylicous.generators.helper_functions import theta
from eddylicous.generators.helper_functions import delta_star
from eddylicous.readers.foamfile_readers import read_points_from_foamfile
from eddylicous.writers.tvmfv_writers import write_points_to_tvmfv
from eddylicous.writers.hdf5_writers import write_points_to_hdf5
from eddylicous.generators.lund_rescaling import lund_generate
from eddylicous.generators.lund_rescaling import lund_rescale_mean_velocity


# Define the command-line arguments
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

readPath = configDict["readPath"]
inflowReadPath = configDict["inflowReadPath"]
writePath = configDict["writePath"]

sampleSurfaceName = configDict["sampleSurfaceName"]
inletPatchName = configDict["inletPatchName"]

reader = configDict["reader"]
inflowReader = configDict["inflowReader"]
writer = configDict["writer"]

hdf5FileName = configDict["hdf5FileName"]

if (reader == "foamFile"):
    dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
else:
    print "ERROR in runLundRescaling.py: unknown reader ", configDict["reader"]
    exit()


# Grab the existing times and sort
if (reader == "foamFile"):
    times = os.listdir(dataDir)
    times = np.sort(times)
else:
    print "ERROR in runLundRescaling.py: unknown reader ", configDict["reader"]
    exit()

# Get the mean profile
if (reader == "foamFile"):
    uMean = np.append(np.zeros((1, 1)),
                      np.genfromtxt(os.path.join(readPath, "postProcessing",
                                                 "collapsedFields", "240",
                                                 "UMean_X.xy"))[:, 1])
else:
    print "ERROR in runLundRescaling.py: unknown reader ", configDict["reader"]
    exit()

nPointsY = uMean.size

# Read grid for the recycling plane
if (reader == "foamFile"):
    [pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
        os.path.join(dataDir, times[0], sampleSurfaceName, "faceCentres"),
        nPointsY=nPointsY)
else:
    print "ERROR in runLundRescaling.py: unknown reader ", configDict["reader"]
    exit()

[nPointsY, nPointsZ] = pointsY.shape

# Read grid for the inflow plane
if (inflowReader == "foamFile"):
    [pointsYInfl, pointsZInfl, yIndInfl, zIndInfl] = read_points_from_foamfile(
        os.path.join(inflowReadPath, inletPatchName, "faceCentres"),
        addZeros=False)
else:
    print "ERROR in runLundRescaling.py: unknown reader ", \
          configDict["inflowReader"]
    exit()

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
uTauInfl = Ue*np.sqrt(cfInfl/2)
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
t = t0

# Get the grid points along y as 1d arrays for convenience
yPrec = pointsY[:, 0]
yInfl = pointsYInfl[:, 0] - yOrigin

deltaPrec = delta_99(yPrec, uMean)
ReTauPrec = uTauPrec*deltaPrec/nuPrec

# Outer scale coordinates
etaPrec = yPrec/deltaPrec
etaInfl = yInfl/deltaInfl

# Inner scale coordinates
yPlusPrec = yPrec*uTauPrec/nuPrec
yPlusInfl = yInfl*uTauInfl/nuInfl

# Points containing the boundary layer at the inflow plane
nInfl = 0
for i in xrange(etaInfl.size):
    if etaInfl[i] <= etaPrec[-1]:
        nInfl += 1

# Points where inner scaling will be used
nInner = 0
for i in xrange(etaInfl.size):
    if etaInfl[i] <= 0.7:
        nInner += 1


# Simple sanity checks
if (deltaInfl > yInfl[-1]):
    print "ERROR in runLundRescaling.py. Desired delta_99 is \
          larger then maximum y."
    exit()

if ReTauInfl > ReTauPrec:
    print "WARNING: Re_tau in the precursor is lower than in the desired TBL"


# Write points and modify writePath appropriatly
if (writer == "tvmfv"):
    write_points_to_tvmfv(os.path.join(writePath, "constant", "boundaryData",
                                       "points"),
                          pointsYInfl, pointsZInfl, xOrigin)
    writePath = os.path.join(writePath, "constant", "boundaryData",
                             inletPatchName)
elif (writer == "hdf5"):
    writePath = os.path.join(writePath, hdf5FileName)
    write_points_to_hdf5(writePath, pointsYInfl, pointsZInfl, xOrigin)
else:
    print "ERROR in runLundRescaling.py. Unknown writer ", configDict["writer"]
    exit()

uMeanInfl = lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMean,
                                       nInfl, nInner,
                                       etaInfl, yPlusInfl, nPointsZInfl,
                                       Ue, U0, gamma)

ReThetaInfl = theta(yInfl, uMeanInfl[:, 0])*Ue/nuInfl
ReDeltaStarInfl = delta_star(yInfl, uMeanInfl[:, 0])*Ue/nuInfl

# Generate the inflow fields

print "Generating the inflow fields."
lund_generate(reader, dataDir,
              writer, writePath,
              dt, t0,
              uMean, uMeanInfl,
              etaPrec, yPlusPrec, pointsZ,
              etaInfl, yPlusInfl, pointsZInfl,
              nInfl, nInner, gamma,
              yInd, zInd,
              surfaceName=sampleSurfaceName, times=times)
print "Done."
