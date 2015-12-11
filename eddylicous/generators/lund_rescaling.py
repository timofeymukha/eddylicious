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

    This function rescales the mean velocity profile taken from
    the precursor simulation using Lund et al's rescaling.


    Parameters
    ----------
        etaPrec : 1d ndarray
            The values of eta for the corresponding values
            of the mean velocity from the precursor.
        yPlusPrec : ndarray
            The values of y+ for the corresponding values
            of the mean velocity from the precursor.
        uMeanPrec : 1d ndarray
            The values of the mean velocity from the precursor.
        nInfl : int
            The amount of points in the wall-normal direction
            that contain the boundary layer at the inflow
            boundary. That is, for points beyound nInfl, Ue will
            be prescribed.
        nInner : int
            The amount of points where inner rescaling should be
            considered. For points beyound nInner, the outer
            rescaling only we be computed. The relaxes the demand
            on Re_tau for the precursor.
        etaInfl : 1d ndarray
            The values of eta for the mesh points at the inflow
            boundary.
        yPlusInfl : 1d ndarray
            The values of y+ for the meshpoints at the inflow
            boundary.
        nPointsZInfl : int
            The amount of points in the spanwise direction for
            the inflow boundary.
        Ue : float
            The freestream velocity.
        U0 : float
            The centerline velocity for the precursor.
        gamma : float
            The ration of the friction velocities in the inflow
            boundary layer and the precursor.


    Returns
    -------
        ndarray
            A 2d ndarray with the values of the mean velocity.
            As expected, the values only vary in the y direction.
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
    """Rescale the fluctuations of velocity using Lund et al's
    rescaling.

    This function rescales the fluctuations of the three
    components of the velocity field taken from the precursor
    simulation using Lund et al's rescaling.


    Parameters
    ----------
        etaPrec : ndarray
            The values of eta for the corresponding values
            of the mean velocity from the precursor.
        yPlusPrec : ndarray
            The values of y+ for the corresponding values
            of the mean velocity from the precursor.
        pointsZ : ndarray
            A 2d array containing the values of z for the points
            of the precuror mesh.
        uPrimeX : ndarray
            A 2d array containing the values of the fluctuations
            of the x component of velocity.
        uPrimeY : ndarray
            A 2d array containing the values of the fluctuations
            of the y component of velocity.
        uPrimeZ : ndarray
            A 2d array containing the values of the fluctuations
            of the z component of velocity.
        gamma : float
            The ration of the friction velocities in the inflow
            boundary layer and the precursor.
        etaInfl : 1d ndarray
            The values of eta for the mesh points at the inflow
            boundary.
        yPlusInfl : 1d ndarray
            The values of y+ for the meshpoints at the inflow
            boundary.
        pointsZInfl : int
            A 2d array containing the values of z for the points
            of the inflow boundary.
        nInfl : int
            The amount of points in the wall-normal direction
            that contain the boundary layer at the inflow
            boundary. That is, for points beyound nInfl, 0 will
            be prescribed.
        nInner : int
            The amount of points where inner rescaling should be
            considered. For points beyound nInner, the outer
            rescaling only we be computed. The relaxes the demand
            on Re_tau for the precursor.


    Returns
    -------
        List of ndarrays
            The list contains three items, each a 2d ndarray.
            The first array contains the rescaled fluctuations of
            the x component of veloicty. The second -- of the y
            component of velocity. The third -- of the z component
            of velocity.
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


def lund_generate(reader, readPath,
                  writer, writePath,
                  dt, t0,
                  uMeanPrec, uMeanInfl,
                  etaPrec, yPlusPrec, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  nInfl, nInner, gamma,
                  yInd, zInd,
                  surfaceName="None", times="None"):
    """Generate the files with the inflow velocity using Lund's rescaling.

        This function will use Lund et al's rescaling in order to
        generate velocity fields for the inflow boundary.
        The rescaling for the mean profile should be done before-
        hand and is one of the input parameters for this funcion.


        Parameters
        ----------
        reader : str
            The name of the reader which will be used to read the
            velocity field from the precursor simulation.
        readPath : str
            The path for the reader.
        writer: str
            The writer that will be used to save the values of the
            velocity field.
        writePath : str
            The path for the writer
        dt : float
            The time-step to be used in the simulation. This will be
            used to associate a time-value with the produced velocity
            fields.
        t0 : float
            The starting time to be used in the simulation. This will
            be used to associate a time-value with the produced velocity
            fields.
        uMeanPrec : 1d ndarray
            The values of the mean velocity from the precursor.
        uMeanInfl : 1d ndarray
            The values of the mean velocity for the inflow boundary
            layer.
        etaPrec : ndarray
            The values of eta for the corresponding values
            of the mean velocity from the precursor.
        yPlusPrec : ndarray
            The values of y+ for the corresponding values
            of the mean velocity from the precursor.
        pointsZ : ndarray
            A 2d array containing the values of z for the points
            of the precuror mesh.
        etaInfl : 1d ndarray
            The values of eta for the mesh points at the inflow
            boundary.
        yPlusInfl : 1d ndarray
            The values of y+ for the meshpoints at the inflow
            boundary.
        pointsZInfl : int
            A 2d array containing the values of z for the points
            of the inflow boundary.
        nInfl : int
            The amount of points in the wall-normal direction
            that contain the boundary layer at the inflow
            boundary. That is, for points beyound nInfl, 0 will
            be prescribed.
        nInner : int
            The amount of points where inner rescaling should be
            considered. For points beyound nInner, the outer
            rescaling only we be computed. The relaxes the demand
            on Re_tau for the precursor.
        gamma : float
            The ration of the friction velocities in the inflow
            boundary layer and the precursor.
        yInd : ndarray
            The sort indices for sorting the read velocity field.
            This is needed when some sorting is performed when the
            mesh points are read and turned into ordered 2d arrays.
            Them the exact same sorting should be applyed to the
            velocity fields.
        zInd : ndarray
            Same as yInd, but for the sorting of the z values.
        surfaceName : str, optional
            For the foamFile reader, the name of the surface used for
            sampling the velocity values.
        times : list of floats, optional
            For the foamFile reader, the times for which the velocity
            field was sampled in the precursor simulation.
    """

    t = t0
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
