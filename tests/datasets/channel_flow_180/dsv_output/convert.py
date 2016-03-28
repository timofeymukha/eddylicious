import numpy as np
import os

"""
A small script to convert dsv as taken from foam-files
into 2d NumPy arrays
"""

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

times = [ name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name)) ]
times = np.sort(times)

for i in xrange(len(times)):
    u = np.genfromtxt(os.path.join(times[i], "U"))
    # Sort along y
    u[:, 0] = u[yInd, 0]
    u[:, 1] = u[yInd, 1]
    u[:, 2] = u[yInd, 2]

    # Reshape to 2d
    uX = np.copy(np.reshape(u[:, 0], (-1, nPointsZ)))
    uY = np.copy(np.reshape(u[:, 1], (-1, nPointsZ)))
    uZ = np.copy(np.reshape(u[:, 2], (-1, nPointsZ)))

    # Sort along z
    for j in xrange(uX.shape[0]):
        uX[j, :] = uX[j, zInd[j, :]]
        uY[j, :] = uY[j, zInd[j, :]]
        uZ[j, :] = uZ[j, zInd[j, :]]
    np.save(os.path.join(times[i], "uX"), uX)    
    np.save(os.path.join(times[i], "uY"), uY)    
    np.save(os.path.join(times[i], "uZ"), uZ)    
