import os
import numpy as np
from sys import exit
from mpi4py import MPI 
from scipy.interpolate import interp1d
from .helper_functions import chunks_and_offsets
from .helper_functions import chauhan_U_composite
from .helper_functions import epsilon_ReT
from .lund_rescaling import lund_rescale_fluctuations
from eddylicious.readers.foamfile_readers import read_u_from_foamfile
from eddylicious.writers.tvmfv_writers import write_u_to_tvmfv
from eddylicious.writers.hdf5_writers import write_u_to_hdf5

__all__ = ["mod_lund_rescale_mean_velocity", "mod_lund_generate"]

""" Our new method for generating the inflow profile

"""


def mod_lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMeanPrec,
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
    # Chauhan profile for the TBL
    lawTbl = interp1d(etaInfl[:nInfl], 
                      chauhan_U_composite(yPlusInfl[:nInfl],
                                          etaInfl[:nInfl], 0.26))

    
    # Reynolds-T profile for channel flow
    ReTauPrec = yPlusPrec[-1]
    c0Kappa = 0.4326
    c1Kappa = 4.2321e+3
    mKappa = 1.9228
    kappaPrec = c0Kappa + c1Kappa*pow(ReTauPrec, -mKappa) 

    c0APlus = 28.1539
    c1APlus = 2.5345e+6
    mAPlus = 2.2619
    APlusPrec = c0APlus + c1APlus*pow(ReTauPrec, -mAPlus) 

    epsilon = epsilon_ReT(etaPrec, ReTauPrec, kappaPrec, APlusPrec)

    lawChannel = np.zeros(etaPrec.shape) 

    for i in xrange(lawChannel.shape[0]):
         lawChannel[i] = np.trapz((1-etaPrec[:i+1])/(1+epsilon[:i+1]),
                                  x=etaPrec[:i+1])
    lawChannel *= 550
    lawChannel = interp1d(etaPrec, lawChannel)


    uMeanInterp = interp1d(etaPrec, uMeanPrec)

    uMeanInfl = np.zeros(etaInfl.shape)
    uMeanInfl[:nInfl] = gamma*lawTbl(etaInfl[:nInfl])/lawChannel(etaInfl[:nInfl])
    uMeanInfl[:nInfl] *= uMeanInterp(etaInfl[:nInfl])

    uMeanInfl[nInfl:] = uMeanInfl[nInfl-1]
    uMeanInfl = np.ones((etaInfl.size, nPointsZInfl))*uMeanInfl[:, np.newaxis]
    return uMeanInfl


def mod_lund_generate(reader, readPath,
                      writer, writePath,
                      dt, t0, tEnd, timePrecision,
                      uMeanPrec, uMeanInfl,
                      etaPrec, yPlusPrec, pointsZ,
                      etaInfl, yPlusInfl, pointsZInfl,
                      nInfl, nInner, gamma,
                      yInd, zInd,
                      surfaceName="None", times="None"):
    """To be documented.

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
           The path for the writer.
       dt : float
           The time-step to be used in the simulation. This will be
           used to associate a time-value with the produced velocity
           fields.
       t0 : float
           The starting time to be used in the simulation. This will
           be used to associate a time-value with the produced velocity
       timePrecision : int
           Number of points after the decimal to keep for the time value.
       tEnd : float
           The ending time for the simulation.
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

    # Grab info regarding parallelization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()

    # Get the total amount of rescalings to be done
    size = int((tEnd-t0)/dt+1)
    
    # Calculate the amount of rescaligns each processor is responsible for
    [chunks, offsets] = chunks_and_offsets(nProcs, size)

    # Perform the rescaling
    for i in xrange(chunks[rank]):
        t = t0 + dt*i + dt*int(offsets[rank])
        t = float(("{0:."+str(timePrecision)+"f}").format(t))
        position = int(offsets[rank]) + i

        if (rank == 0):
            print "     Rescaled about", i/float(chunks[rank])*100, "%" 

        # Read U data
        flag = False 
        if (reader == "foamFile"):
            if (position < len(times)):
                readPosition = position
            else:
                if (not flag):
                    print "Warning: precursor database small than required \
                    number of time-steps"
                    flag = True
                readPosition = position - np.floor(position/len(times))*len(times)
                
            [U_X, U_Y, U_Z] = read_u_from_foamfile(
                os.path.join(readPath, times[readPosition], surfaceName,
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

        # Write
        if (writer == "tvmfv"):
            write_u_to_tvmfv(writePath, t, UInfl)
        elif (writer == "hdf5"):
            write_u_to_hdf5(writePath, t, UInfl, position, size)
        else:
            print "ERROR in mod_lund_generate(). Unknown writer ", writer
            exit()

