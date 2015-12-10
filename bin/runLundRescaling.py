import os
import shutil
import numpy as np
from sys import exit
import matplotlib.pyplot as plt


#channelCaseDir = "../../../../project-2.3.1/run/channel_flow/channel550"
channelCaseDir = "../channelPronkRefined"
mainCaseDir = "../../ramp/bll/bll_coarse"

surfaceName = "inletSurface"
mainPatchName = "inletBot"

dataDir = os.path.join(channelCaseDir, "postProcessing", "sampledSurface")
#dataDir = path.join("..", "databases", "channel395_1e-1")
boundaryDataDir = os.path.join(mainCaseDir, "constant", "boundaryData", mainPatchName)

# Grab the existing times and sort
times = os.listdir(dataDir)

# Get the mean profile
UMean = np.append( np.zeros((1,1)), np.genfromtxt(os.path.join(channelCaseDir,
    "postProcessing", "collapsedFields", "240", "UMean_X.xy"))[:,1])
UMeanTbl = np.genfromtxt("../tblDNS/Schlatter/vel_670_dns.prof")[:,0:3]


NPointsY = UMean.size

# Read grid for the recycling plane 
[pointsY, pointsZ, yInd, zInd ]= read_points_from_foamfile(
    os.path.join(dataDir, times[0], surfaceName, "faceCentres"), 
    NPointsY=NPointsY)

[NPointsY, NPointsZ] = pointsY.shape

# Read grid for the inflow plane 
[pointsYInfl, pointsZInfl, yIndInfl, zIndInfl] = read_points_from_foamfile(
    os.path.join(mainCaseDir, "postProcessing", "surfaces", "0", mainPatchName, "faceCentres"),
    #path.join(dataDir, times[0], surfaceName, "faceCentres"),
    addZeros=False)
    
[NPointsYInfl, NPointsZInfl] = pointsYInfl.shape


# TBL AND CHANNEL FLOW PARAMETERS 
nuInfl = 0.0000023 
nuRecy = 0.001937 
 

# Boundary layer thickneses
deltaInfl = 0.85*31.5e-3 
deltaRecy = 1

# Freestream velocity, and centerline velocity
U0 = np.max(UMean) 
Ue = 1 

# Reynolds number based on delta
ReDelta = Ue*deltaInfl/nuInfl

# Cf at the inflow, to get u_tau
CfInfl = 0.02*pow(1.0/ReDelta, 1.0/6)

# u_tau 
uTauInfl = Ue*np.sqrt(CfInfl/2)
uTauRecy = np.sqrt(nuRecy*UMean[1]/pointsY[1,0])

ReTauInfl = uTauInfl*deltaInfl/nuInfl
ReTauRecy = uTauRecy*deltaRecy/nuRecy

# gamma, the ratio of friction velocities
gamma = uTauInfl/uTauRecy

# Shift in coordinates
xShift = -0.23121
yShift = 0.127

# Get the grid points along y as 1d arrays for convenience 
yRecy = pointsY[:,0]
yInfl = pointsYInfl[:,0]-yShift

deltaRecy = delta_99(yRecy, UMean)
#deltaInfl = deltaRecy

# Outer scale coordinates
etaRecy = yRecy/deltaRecy 
etaInfl = yInfl/deltaInfl

# Inner scale coordinates
yPlusRecy = yRecy*uTauRecy/nuRecy
yPlusInfl = yInfl*uTauInfl/nuInfl


# Sanity checks
if (deltaInfl > yInfl[-1]):
    print "FATAL ERROR: desired delta_99 is larger then the upper boundary of the mesh"
    exit()

if ReTauInfl > ReTauRecy:
    print "WARNING: Re_tau in the channel flow is lower than in the desired TBL\n"

# Points containing the boundary layer at the inflow plane

    # Headers to add to the files
scalarHeader = "FoamFile\n{\nversion 2.0;\nformat ascii;\nclass scalarAverageField;\nobject values;\n}\n\n0\n"


write_points_to_tvmfv(os.path.join(boundaryDataDir, "points"), pointsYInfl, pointsZInfl, xShift)

UMeanInfl = lund_rescale_mean_velocity(etaRecy, yPlusRecy, UMean,
                                        etaInfl, yPlusInfl)


ReThetaInfl = theta(yInfl, UMeanInfl[:,0])*Ue/nuInfl

#exit()

dt = 5e-4 
t = 0
timePrecision = 5

times = np.sort(times)


#lund_generate(read_u_from_foamfile, dataDir, 
#             write_u_to_tvmfv, boundaryDataDir, 
#             times, dt,
#             UMean, UMeanInfl,
#             etaRecy, yPlusRecy, pointsZ,
#             etaInfl, yPlusInfl, pointsZInfl)


WVals = W(etaInfl)[:,np.newaxis]

N = 45000

print "Loading database..."
UPrime_X = np.load("dbtest/UPrime_X.npy")[:,:,:N]
UPrime_Y = np.load("dbtest/UPrime_Y.npy")[:,:,:N]
UPrime_Z = np.load("dbtest/UPrime_Z.npy")[:,:,:N]

#interp_X = [interp2d(pointsZ[0,:]/pointsZ[0,-1], etaRecy, UPrime_X[:,:,i])\
#        for i in xrange(N)]
#plusInterp_X = [interp2d(pointsZ[0,:]/pointsZ[0,-1], yPlusRecy, UPrime_X[:,:,i])\
#        for i in xrange(N)]

#UPrime_X_infl = [(1-WVals[:NInfl,:])*gamma*plusInterp_X[i](pointsZInfl[0,:]/pointsZInfl[0,-1], 
#    yPlusInfl[0:NInfl])+ WVals[:NInfl,:]*gamma*interp_X[i](pointsZInfl[0,:]/pointsZInfl[0,-1], 
#    etaInfl[0:NInfl]) for i in xrange(len(plusInterp_X))]

#UPrime_X_infl = np.dstack(UPrime_X_infl)

#UPrime2Mean_XX  = np.mean(np.var(UPrime_X, axis=2), axis=1)
#UPrimeMean_X  = np.mean(np.mean(UPrime_X, axis=2), axis=1)

#UPrime2Mean_XX_infl  = np.mean(np.var(UPrime_X_infl, axis=2), axis=1)
#UPrimeMean_X_infl  = np.mean(np.mean(UPrime_X_infl, axis=2), axis=1)

#del interp_X
#del plusInterp_X
#del UPrime_X
#gc.collect()

print "Rescaling!"


printCounter = -1
for timeI in xrange(len(times)):
    printCounter += 1
    if (printCounter == 100):
        print timeI
        printCounter = 0

    if timeI > 40000:
        break

    [UPrime_X_infl, UPrime_Y_infl, UPrime_Z_infl] = \
                lund_rescale_fluctuations(etaRecy, yPlusRecy, pointsZ,
                                         UPrime_X[:,:,timeI],
                                         UPrime_Y[:,:,timeI], 
                                         UPrime_Z[:,:,timeI],
                                         etaInfl, yPlusInfl, pointsZInfl)

    # Combine and flatten
    U_infl_X = np.reshape(UPrime_X_infl+UMeanInfl, (UPrime_X_infl.size, -1), order='F')
    U_infl_Y = np.reshape(UPrime_Y_infl, (UPrime_X_infl.size, -1), order='F')
    U_infl_Z = np.reshape(UPrime_Z_infl, (UPrime_X_infl.size, -1), order='F')

    UInfl = np.concatenate((U_infl_X, U_infl_Y, U_infl_Z), axis = 1)

    write_u_to_tvmfv(boundaryDataDir, t, UInfl)

    t += dt

