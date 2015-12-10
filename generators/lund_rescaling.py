import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d


def lund_rescale_mean_velocity(etaRecy, yPlusRecy, uMeanRecy,
                               etaInfl, yPlusInfl, nPointsZInfl,
                               Ue, U0, gamma):
    """Rescale the mean velocity profile using Lunds rescaling.

    Returns a 2d numpy array with the values of mean velocity.
    """

# Points containing the boundary layer at the inflow plane
    nInfl = 0
    for i in xrange(etaInfl.size):
        if etaInfl[i] <= etaRecy[-1]:
            nInfl += 1

# Points where inner scaling will be used
    nInner = 0
    for i in xrange(etaInfl.size):
        if etaInfl[i] <= 0.7:
            nInner += 1

    uMeanInterp = interp1d(etaRecy, uMeanRecy)
    uMeanInterpPlus = interp1d(yPlusRecy, uMeanRecy)

    uMeanInner = np.append(gamma*uMeanInterpPlus(yPlusInfl[0:nInner]),
                           np.zeros(nInfl-nInner))
    uMeanOuter = gamma*uMeanInterp(etaInfl[0:nInfl]) + Ue - gamma*U0

    uMeanInfl = np.zeros(etaInfl.shape)
    uMeanInfl[0:nInfl] = uMeanInner*(1-blending_function(etaInfl[0:nInfl])) + \
        uMeanOuter*blending_function(etaInfl[0:nInfl])
    uMeanInfl[nInfl:] = Ue
    uMeanInfl = np.ones((etaInfl.size, nPointsZInfl))*uMeanInfl[:, np.newaxis]
    return uMeanInfl


def lund_rescale_fluctuations(etaRecy, yPlusRecy, pointsZ,
                              uPrimeX, uPrimeY, uPrimeZ, gamma,
                              etaInfl, yPlusInfl, pointsZInfl, nInfl):
    """Rescale the fluctuations using Lund's rescaling.

    Returns a list with 3 items: numpy arrays for each of
    the components of the fluctuations.
    """

    uPrimeXInfl = np.zeros(pointsZInfl.shape)
    uPrimeYInfl = np.zeros(pointsZInfl.shape)
    uPrimeZInfl = np.zeros(pointsZInfl.shape)

    uPrimeXInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaRecy, uPrimeX)
    uPrimeYInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaRecy, uPrimeY)
    uPrimeZInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaRecy, uPrimeZ)

    uPrimeXPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusRecy,
                                 uPrimeX)
    uPrimeYPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusRecy,
                                 uPrimeY)
    uPrimeZPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusRecy,
                                 uPrimeZ)

    uPrimeXInner = \
        gamma*uPrimeXPlusInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                yPlusInfl[0:nInfl])
    uPrimeYInner = \
        gamma*uPrimeYPlusInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                yPlusInfl[0:nInfl])
    uPrimeZInner = \
        gamma*uPrimeZPlusInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                yPlusInfl[0:nInfl])

    uPrimeXOuter = gamma*uPrimeXInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                       etaInfl[0:nInfl])
    uPrimeYOuter = gamma*uPrimeYInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                       etaInfl[0:nInfl])
    uPrimeZOuter = gamma*uPrimeZInterp(pointsZInfl[0, :]/pointsZInfl[0, -1],
                                       etaInfl[0:nInfl])

    uPrimeXInfl[0:nInfl] = \
        uPrimeXInner*(1-blending_function(etaInfl[0:nInfl]))[:, np.newaxis] + \
        uPrimeXOuter*blending_function(etaInfl[0:nInfl])[:, np.newaxis]
    uPrimeYInfl[0:nInfl] = \
        uPrimeYInner*(1-blending_function(etaInfl[0:nInfl]))[:, np.newaxis] + \
        uPrimeYOuter*blending_function(etaInfl[0:nInfl])[:, np.newaxis]
    uPrimeZInfl[0:nInfl] = \
        uPrimeZInner*(1-blending_function(etaInfl[0:nInfl]))[:, np.newaxis] + \
        uPrimeZOuter*blending_function(etaInfl[0:nInfl])[:, np.newaxis]

    return [uPrimeXInfl, uPrimeYInfl, uPrimeZInfl]


def lund_generate(reader, readPath, surfaceName,
                  writer, writePath,
                  times, dt,
                  uMeanRecy, uMeanInfl,
                  etaRecy, yPlusRecy, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  yInd, zInd):
    """Generate the files with the inflow velocity using
    Lund's rescaling.
    """

    t = 0
    for timeI in xrange(len(times)):
        print timeI

        # Read U data
        [U_X, U_Y, U_Z] = reader(
            os.path.join(readPath, times[timeI], surfaceName,
                         "vectorField", "U"),
            uMeanRecy.size,
            yInd,
            zInd)

        # Claculate UPrime
        uPrimeX = U_X - uMeanRecy[:, np.newaxis]
        uPrimeY = U_Y
        uPrimeZ = U_Z

        [uPrimeXInfl, uPrimeYInfl, uPrimeZInfl] = \
            lund_rescale_fluctuations(
                etaRecy, yPlusRecy, pointsZ,
                uPrimeX, uPrimeY, uPrimeZ,
                etaInfl, yPlusInfl, pointsZInfl)

        # Combine and flatten
        UInfl_X = np.reshape(uPrimeXInfl+uMeanInfl, (uPrimeXInfl.size, -1),
                             order='F')
        UInfl_Y = np.reshape(uPrimeYInfl, (uPrimeXInfl.size, -1), order='F')
        UInfl_Z = np.reshape(uPrimeZInfl, (uPrimeXInfl.size, -1), order='F')

        UInfl = np.concatenate((UInfl_X, UInfl_Y, UInfl_Z), axis=1)

        writer(writePath, t, UInfl)

        t += dt
