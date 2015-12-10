import os
import numpy as np
from sys import exit
from eddylicous.generators.helper_functions import delta_99
from eddylicous.generators.helper_functions import theta
from eddylicous.generators.helper_functions import delta_star
from eddylicous.readers.foamfile_readers import read_points_from_foamfile
from eddylicous.readers.foamfile_readers import read_u_from_foamfile
from eddylicous.writers.tvmfv_writers import write_u_to_tvmfv
from eddylicous.writers.tvmfv_writers import write_points_to_tvmfv
from eddylicous.generators.lund_rescaling import lund_generate
from eddylicous.generators.lund_rescaling import lund_rescale_mean_velocity


# Read the config file into a dictionary
configDict = {}

configFile = open('testConfig', mode='r')

for line in configFile:
    if (line[0] == '#') or (line == '\n'):
        continue
    configDict[line.split()[0]] = line.split()[1]

readPath = configDict["readPath"]
mainCaseDir = configDict["writePath"]

sampleSurfaceName = configDict["samleSurfaceName"]
inletPatchName = configDict["inletPatchName"]

reader = configDict["reader"]
inflowReader = configDict["inflowReader"]
writer = configDict["writer"]

if (reader == "foamFile"):
    dataDir = os.path.join(readPath, "postProcessing", "sampledSurface")
else:
    print "ERROR: unknown reader ", configDict["reader"]
    exit()

if (writer == "tvmfv"):
    boundaryDataDir = os.path.join(mainCaseDir, "constant",
                                   "boundaryData", inletPatchName)
else:
    print "ERROR. Unknown writer ", configDict["writer"]
    exit()

# Grab the existing times and sort
if (reader == "foamFile"):
    times = os.listdir(dataDir)
    times = np.sort(times)
else:
    print "ERROR: unknown reader ", configDict["reader"]
    exit()

# Get the mean profile
if (reader == "foamFile"):
    uMean = np.append(np.zeros((1, 1)),
                      np.genfromtxt(os.path.join(readPath, "postProcessing",
                                                 "collapsedFields", "240",
                                                 "UMean_X.xy"))[:, 1])
else:
    print "ERROR: unknown reader ", configDict["reader"]
    exit()

nPointsY = uMean.size

# Read grid for the recycling plane
if (reader == "foamFile"):
    [pointsY, pointsZ, yInd, zInd] = read_points_from_foamfile(
        os.path.join(dataDir, times[0], sampleSurfaceName, "faceCentres"),
        nPointsY=nPointsY)
else:
    print "ERROR: unknown reader ", configDict["reader"]
    exit()

[nPointsY, nPointsZ] = pointsY.shape

# Read grid for the inflow plane
if (inflowReader == "foamFile"):
    [pointsYInfl, pointsZInfl, yIndInfl, zIndInfl] = read_points_from_foamfile(
        os.path.join(mainCaseDir, "postProcessing", "surfaces", "0",
                     inletPatchName, "faceCentres"),
        addZeros=False)
else:
    print "ERROR: unknown reader ", configDict["inflowReader"]
    exit()

[nPointsYInfl, nPointsZInfl] = pointsYInfl.shape

# Define and compute the parameters of the boundary layer

# Viscosity
nuInfl = float(configDict["nuInflow"])
nuPrec = float(configDict["nuPrecursor"])


# Boundary layer thickneses
deltaInfl = float(configDict["delta"])

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
xOrigen = float(configDict["xOrigen"])
yOrigen = float(configDict["yOrigen"])

# Time-step and initial time for the writer
dt = float(configDict["dt"])
t0 = float(configDict["t0"])
t = t0

# Get the grid points along y as 1d arrays for convenience
yPrec = pointsY[:, 0]
yInfl = pointsYInfl[:, 0] - yOrigen

deltaPrec = delta_99(yPrec, uMean)
ReTauPrec = uTauPrec*deltaPrec/nuPrec

# Outer scale coordinates
etaPrec = yPrec/deltaPrec
etaInfl = yInfl/deltaInfl

# Inner scale coordinates
yPlusPrec = yPrec*uTauPrec/nuPrec
yPlusInfl = yInfl*uTauInfl/nuInfl


# Sanity checks
if (deltaInfl > yInfl[-1]):
    print "ERROR. Desired delta_99 is larger then maximum y."
    exit()

if ReTauInfl > ReTauPrec:
    print "WARNING: Re_tau in the precursor is lower than in the desired TBL"


if (writer == "tvmfv"):
    write_points_to_tvmfv(os.path.join(boundaryDataDir, "points"), pointsYInfl,
                          pointsZInfl, xOrigen)
else:
    print "ERROR. Unknown writer ", configDict["writer"]
    exit()

uMeanInfl = lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMean,
                                       etaInfl, yPlusInfl)

ReThetaInfl = theta(yInfl, uMeanInfl[:, 0])*Ue/nuInfl
ReDeltaStarInfl = delta_star(yInfl, uMeanInfl[:, 0])*Ue/nuInfl

lund_generate(read_u_from_foamfile, dataDir,
              write_u_to_tvmfv, boundaryDataDir,
              times, dt,
              uMean, uMeanInfl,
              etaPrec, yPlusPrec, pointsZ,
              etaInfl, yPlusInfl, pointsZInfl)
