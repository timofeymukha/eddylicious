import os
from os import path
import shutil
import numpy as np
from sys import exit
import matplotlib.pyplot as plt


def read_points_from_foamfile(path, addZeros=1, NPointsY=0, delta=1):
    with file(path) as pointsFile:
        points = [line.rstrip(')\n') for line in pointsFile]            

    points = [line.lstrip('(') for line in points]
    NPoints = int(points[1])
    points = points[3:-1]
    points = np.genfromtxt(points)[:,1:]

# Sort the points
        #Sort along y first
    yInd = np.argsort(points[:,0])
    points[:,0] = points[yInd,0] 
    points[:,1] = points[yInd,1] 


        # Find the number of points along z
    NPointsZ = 0
    for i in xrange(points[:,0].size):
        if points[i,0] == points[0,0]:
            NPointsZ += 1
        else:
            break

        # Reshape into a 2d array
    pointsY = np.copy(np.reshape(points[:,0], (-1,NPointsZ)))
    pointsZ = np.copy(np.reshape(points[:,1], (-1,NPointsZ)))

        # For each y order the points in z

    zInd = np.zeros(pointsZ.shape, dtype=np.int)

    for i in xrange(pointsZ.shape[0]):
        zInd[i,:] = np.argsort(pointsZ[i,:])
        pointsZ[i,:] = pointsZ[i,zInd[i,:]]


        # Add points at y = 0
    if addZeros:
        pointsY = np.append(np.zeros( (1, pointsY.shape[1]) ), pointsY, axis=0)
        pointsZ = np.append(np.array([pointsZ[0,:]]), pointsZ, axis=0)

        # Cap the points, to include ony one half of the channel
        # Makes y=delta the last point
    if NPointsY:
        pointsY = pointsY[:NPointsY,:]
        pointsZ = pointsZ[:NPointsY,:]

        pointsY[-1,:] = delta 
    else:
        NPointsY = pointsY.shape[0]

    return [pointsY, pointsZ, yInd, zInd]

def read_U_from_foamfile(path, NPointsY, yInd, zInd):

    with file(path) as UFile:
        U = [line.rstrip(')\n') for line in UFile]

    U = [line.lstrip('(') for line in U]
    U = U[3:-1]
    U = np.genfromtxt(U)

    # Sort along y
    U[:,0] = U[yInd,0]
    U[:,1] = U[yInd,1]
    U[:,2] = U[yInd,2]

    # Reshape to 2d
    U_X = np.copy(np.reshape(U[:,0], (-1,NPointsZ)))
    U_Y = np.copy(np.reshape(U[:,1], (-1,NPointsZ)))
    U_Z = np.copy(np.reshape(U[:,2], (-1,NPointsZ)))

    # Sort along z
    for i in xrange(U_X.shape[0]):
        U_X[i,:] = U_X[i,zInd[i,:]]
        U_Y[i,:] = U_Y[i,zInd[i,:]]
        U_Z[i,:] = U_Z[i,zInd[i,:]]
    
    # Prepend with zeros
    U_X = np.append(np.zeros( (1, pointsY.shape[1]) ), U_X, axis=0)
    U_Y = np.append(np.zeros( (1, pointsY.shape[1]) ), U_Y, axis=0)
    U_Z = np.append(np.zeros( (1, pointsY.shape[1]) ), U_Z, axis=0)

    # Interpolate to get data at y=delta
    U_X[NPointsY-1,:] = 0.5*(U_X[NPointsY-2,:] + U_X[NPointsY,:])
    U_Y[NPointsY-1,:] = 0.5*(U_Y[NPointsY-2,:] + U_Y[NPointsY,:])
    U_Z[NPointsY-1,:] = 0.5*(U_Z[NPointsY-2,:] + U_Z[NPointsY,:])

    # Remove data above y=delta
    U_X = U_X[:NPointsY,:]
    U_Y = U_Y[:NPointsY,:]
    U_Z = U_Z[:NPointsY,:]

    return [U_X, U_Y, U_Z]



channelCaseDir = "../channelPronkRefined"

surfaceName = "inletSurface"

dataDir = path.join(channelCaseDir, "postProcessing", "sampledSurface")
#dataDir = path.join("..", "databases", "channel395_1e-1")

# Grab the existing times and sort
times = os.listdir(dataDir)

# Get the mean profile
uMean = np.append( np.zeros((1,1)), np.genfromtxt(path.join(channelCaseDir, "postProcessing", "collapsedFields", "240", "UMean_X.xy"))[:,1])

NPointsY = uMean.size

# Read grid for the recycling plane 
[pointsY, pointsZ, yInd, zInd ]= read_points_from_foamfile(
    path.join(dataDir, times[0], surfaceName, "faceCentres"), 
    NPointsY=NPointsY)

[NPointsY, NPointsZ] = pointsY.shape

uPrimeX = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeY = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))
uPrimeZ = np.zeros((pointsY.shape[0], pointsY.shape[1], len(times)))

times = np.sort(times)

for timeI in xrange(len(times)):

    #if timeI > 30000:
    #    break

    # Read U data
    [U_X, U_Y, U_Z] = read_U_from_foamfile(
        path.join(dataDir, times[timeI], surfaceName, "vectorField", "U"),
        NPointsY,
        yInd,
        zInd)

    uPrimeX[:,:,timeI] = U_X - uMean[:,np.newaxis]
    uPrimeY[:,:,timeI] = U_Y 
    uPrimeZ[:,:,timeI] = U_Z 

    print timeI

uPrimeX = uPrimeX[:,:,1:]
uPrimeY = uPrimeY[:,:,1:]
uPrimeZ = uPrimeZ[:,:,1:]

np.save("dbtest/uPrimeX", uPrimeX)
np.save("dbtest/uPrimeY", uPrimeY)
np.save("dbtest/uPrimeZ", uPrimeZ)
np.save("dbtest/pointsY", pointsY)
np.save("dbtest/pointsZ", pointsZ)
np.save("dbtest/uMean", uMean)


