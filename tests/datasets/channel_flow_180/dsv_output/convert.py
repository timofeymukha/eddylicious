import numpy as np

points = np.genfromtxt("faceCentres")[:, 1:]
yInd = np.argsort(points[:, 0])
points[:, 0] = points[yInd, 0]
points[:, 1] = points[yInd, 1]

nPointsZ = 72
pointsY = np.copy(np.reshape(points[:, 0], (-1, nPointsZ)))
pointsZ = np.copy(np.reshape(points[:, 1], (-1, nPointsZ)))

zInd = np.zeros(pointsZ.shape, dtype=np.int)

for i in xrange(pointsZ.shape[0]):
    zInd[i, :] = np.argsort(pointsZ[i, :])
    pointsZ[i, :] = pointsZ[i, zInd[i, :]]

np.save("yInd", yInd)    
np.save("zInd", zInd)    
np.save("pointsY", pointsY)    
np.save("pointsZ", pointsZ)    
