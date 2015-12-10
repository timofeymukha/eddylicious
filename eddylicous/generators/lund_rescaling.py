import os
import numpy as np
from sys import exit
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from .helper_functions import blending_function
from ..readers.foamfile_readers import read_u_from_foamfile
from ..writers.tvmfv_writers import write_u_to_tvmfv

"""Function for generating inlfow velocity fields using
Lund et al's rescaling, see

Lund T.S., Wu X., Squires K.D. Generation of turbulent inflow
data for spatially-developing boundary layer simulations.
J. Comp. Phys. 1998; 140:233-58.
"""


def lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMeanPrec,
                               nInfl, nInner,
                               etaInfl, yPlusInfl, nPointsZInfl,
                               Ue, U0, gamma):
    """Rescale the mean velocity profile using Lunds rescaling.

    Returns a 2d numpy array with the values of mean velocity.
    """

    uMeanInterp = interp1d(etaPrec, uMeanPrec)
    uMeanInterpPlus = interp1d(yPlusPrec, uMeanPrec)

    uMeanInner = np.append(gamma*uMeanInterpPlus(yPlusInfl[0:nInner]),
                           np.zeros(nInfl-nInner))
    uMeanOuter = gamma*uMeanInterp(etaInfl[0:nInfl]) + Ue - gamma*U0

    uMeanInfl = np.zeros(etaInfl.shape)
    uMeanInfl[0:nInfl] = uMeanInner*(1-blending_function(etaInfl[0:nInfl])) + \
        uMeanOuter*blending_function(etaInfl[0:nInfl])
    uMeanInfl[nInfl:] = Ue
    uMeanInfl = np.ones((etaInfl.size, nPointsZInfl))*uMeanInfl[:, np.newaxis]
    return uMeanInfl


def lund_rescale_fluctuations(etaPrec, yPlusPrec, pointsZ,
                              uPrimeX, uPrimeY, uPrimeZ, gamma,
                              etaInfl, yPlusInfl, pointsZInfl,
                              nInfl, nInner):
    """Rescale the fluctuations using Lund et al's rescaling.

    Returns a list with 3 items: numpy arrays for each of
    the components of the fluctuations.
    """

    uPrimeXInfl = np.zeros(pointsZInfl.shape)
    uPrimeYInfl = np.zeros(pointsZInfl.shape)
    uPrimeZInfl = np.zeros(pointsZInfl.shape)

    uPrimeXInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaPrec, uPrimeX)
    uPrimeYInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaPrec, uPrimeY)
    uPrimeZInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], etaPrec, uPrimeZ)

    uPrimeXPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusPrec,
                                 uPrimeX)
    uPrimeYPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusPrec,
                                 uPrimeY)
    uPrimeZPlusInterp = interp2d(pointsZ[0, :]/pointsZ[0, -1], yPlusPrec,
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
                  uMeanPrec, uMeanInfl,
                  etaPrec, yPlusPrec, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  nInfl, nInner, gamma,
                  yInd, zInd):
    """Generate the files with the inflow velocity using
    Lund's rescaling.
    """

    t = 0
    for timeI in xrange(len(times)):
        print timeI

        # Read U data
        if (reader == "foamFile"):
            [U_X, U_Y, U_Z] = read_u_from_foamfile(
                os.path.join(readPath, times[timeI], surfaceName,
                             "vectorField", "U"),
                pointsZ.shape[0], pointsZ.shape[1],
                yInd, zInd)

        else:
            print "ERROR. Unknown reader ", reader
            exit()

        # Claculate UPrime
        uPrimeX = U_X - uMeanPrec[:, np.newaxis]
        uPrimeY = U_Y
        uPrimeZ = U_Z

        [uPrimeXInfl, uPrimeYInfl, uPrimeZInfl] = \
            lund_rescale_fluctuations(
                etaPrec, yPlusPrec, pointsZ,
                uPrimeX, uPrimeY, uPrimeZ, gamma,
                etaInfl, yPlusInfl, pointsZInfl,
                nInfl, nInner)

        # Combine and flatten
        UInfl_X = np.reshape(uPrimeXInfl+uMeanInfl, (uPrimeXInfl.size, -1),
                             order='F')
        UInfl_Y = np.reshape(uPrimeYInfl, (uPrimeXInfl.size, -1), order='F')
        UInfl_Z = np.reshape(uPrimeZInfl, (uPrimeXInfl.size, -1), order='F')

        UInfl = np.concatenate((UInfl_X, UInfl_Y, UInfl_Z), axis=1)

        if (writer == "tvmfv"):
            write_u_to_tvmfv(writePath, t, UInfl)
        else:
            print "ERROR. Unknown writer ", writer
            exit()

        t += dt
