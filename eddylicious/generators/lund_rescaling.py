"""Function for generating inflow velocity fields using
Lund et al's rescaling, see

Lund T.S., Wu X., Squires K.D. Generation of turbulent inflow
data for spatially-developing boundary layer simulations.
J. Comp. Phys. 1998; 140:233-58.

"""
import os
import numpy as np
from sys import exit
from mpi4py import MPI 
import h5py as h5py
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from .helper_functions import blending_function
from .helper_functions import chunks_and_offsets
from eddylicious.readers.foamfile_readers import read_u_from_foamfile
from eddylicious.writers.tvmfv_writers import write_velocity_to_tvmfv
from eddylicious.writers.hdf5_writers import write_velocity_to_hdf5

__all__ = ["lund_rescale_mean_velocity", "lund_rescale_fluctuations",
           "lund_generate"]


def lund_rescale_mean_velocity(etaPrec, yPlusPrec, uMeanPrec,
                               nInfl, nInner,
                               etaInfl, yPlusInfl, nPointsZInfl,
                               Ue, U0, gamma):
    """Rescale the mean velocity profile using Lunds rescaling.

    This function rescales the mean velocity profile taken from the
    precursor simulation using Lund et al's rescaling.

    Parameters
    ----------
    etaPrec : 1d ndarray
        The values of eta for the corresponding values of the mean
        velocity from the precursor.
    yPlusPrec : ndarray
        The values of y+ for the corresponding values of the mean
        velocity from the precursor.
    uMeanPrec : 1d ndarray
        The values of the mean velocity from the precursor.
    nInfl : int
        The amount of points in the wall-normal direction that contain
        the boundary layer at the inflow boundary. That is, for points
        beyound nInfl, Ue will
        be prescribed.
    nInner : int
        The amount of points where inner rescaling should be considered.
        For points beyound nInner, the outer rescaling only we be
        computed. The relaxes the demand on Re_tau for the precursor.
    etaInfl : 1d ndarray
        The values of eta for the mesh points at the inflow boundary.
    yPlusInfl : 1d ndarray
        The values of y+ for the meshpoints at the inflow boundary.
    nPointsZInfl : int
        The amount of points in the spanwise direction for the inflow
        boundary.
    Ue : float
        The freestream velocity.
    U0 : float
        The centerline velocity for the precursor.
    gamma : float
        The ration of the friction velocities in the inflow boundary
        layer and the precursor.

    Returns
    -------
    ndarray
        A 2d ndarray with the values of the mean velocity. As expected,
        the values only vary in the y direction.

    """
    assert nInfl > 0
    assert nInner > 0
    assert nPointsZInfl > 0
    assert Ue > 0
    assert U0 > 0
    assert gamma > 0
    assert np.all(etaInfl >= 0)
    assert np.all(etaPrec >= 0)
    assert np.all(yPlusInfl >= 0)
    assert np.all(yPlusPrec >= 0)
    assert np.all(uMeanPrec >= 0)
    assert nInner <= nInfl


    uMeanInterp = interp1d(etaPrec, uMeanPrec)
    uMeanInterpPlus = interp1d(yPlusPrec, uMeanPrec)

    uMeanInner = gamma*uMeanInterpPlus(yPlusInfl[:nInner])
    uMeanOuter = gamma*uMeanInterp(etaInfl[:nInfl]) + Ue - gamma*U0

    uMeanInfl = np.zeros(etaInfl.shape)
    uMeanInfl[:nInner] = uMeanInner*(1-blending_function(etaInfl[:nInner])) + \
        uMeanOuter[:nInner]*blending_function(etaInfl[:nInner])
    uMeanInfl[nInner:nInfl] = uMeanOuter[nInner:nInfl]
    uMeanInfl[nInfl:] = Ue
    uMeanInfl = np.ones((etaInfl.size, nPointsZInfl))*uMeanInfl[:, np.newaxis]

    assert np.all(uMeanInfl >= 0)

    return uMeanInfl


def lund_rescale_fluctuations(etaPrec, yPlusPrec, pointsZ,
                              uPrimeX, uPrimeY, uPrimeZ, gamma,
                              etaInfl, yPlusInfl, pointsZInfl,
                              nInfl, nInner):
    """Rescale the fluctuations of velocity using Lund et al's
    rescaling.

    This function rescales the fluctuations of the three components of
    the velocity field taken from the precursor simulation using Lund et
    al's rescaling.

    Parameters
    ----------
    etaPrec : ndarray
        The values of eta for the corresponding values of the mean
        velocity from the precursor.
    yPlusPrec : ndarray
        The values of y+ for the corresponding values of the mean
        velocity from the precursor.
    pointsZ : ndarray
        A 2d array containing the values of z for the points of the
        precuror mesh.
    uPrimeX : ndarray
        A 2d array containing the values of the fluctuations of the x
        component of velocity.
    uPrimeY : ndarray
        A 2d array containing the values of the fluctuations of the y
        component of velocity.
    uPrimeZ : ndarray
        A 2d array containing the values of the fluctuations of the z
        component of velocity.
    gamma : float
        The ration of the friction velocities in the inflow boundary
        layer and the precursor.
    etaInfl : 1d ndarray
        The values of eta for the mesh points at the inflow boundary.
    yPlusInfl : 1d ndarray
        The values of y+ for the meshpoints at the inflow boundary.
    pointsZInfl : int
        A 2d array containing the values of z for the points of the
        inflow boundary.
    nInfl : int
        The amount of points in the wall-normal direction that contain
        the boundary layer at the inflow boundary. That is, for points
        beyound nInfl, 0 will be prescribed.
    nInner : int
        The amount of points where inner rescaling should be considered.
        For points beyound nInner, the outer rescaling only we be
        computed. The relaxes the demand on Re_tau for the precursor.

    Returns
    -------
    List of ndarrays
        The list contains three items, each a 2d ndarray. The first
        array contains the rescaled fluctuations of the x component of
        veloicty. The second -- of the y component of velocity. The
        third -- of the z component of velocity.
    """

    assert np.all(etaPrec >= 0)
    assert np.all(yPlusPrec >= 0)
    assert np.all(etaInfl >= 0)
    assert np.all(yPlusInfl >= 0)
    assert nInfl > 0
    assert gamma > 0

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
                  dt, t0, tEnd, timePrecision,
                  uMeanPrec, uMeanInfl,
                  etaPrec, yPlusPrec, pointsZ,
                  etaInfl, yPlusInfl, pointsZInfl,
                  nInfl, nInner, gamma,
                  yInd, zInd,
                  surfaceName="None", times="None"):
    """Generate the files with the inflow velocity using Lund's
    rescaling, in parallel.

    This function will use Lund et al's rescaling in order to generate
    velocity fields for the inflow boundary. The rescaling for the mean
    profile should be done before- hand and is one of the input
    parameters for this function.

    Parameters
    ----------
    reader : str
        The name of the reader which will be used to read the velocity
        field from the precursor simulation.
    readPath : str
        The path for the reader.
    writer: str
        The writer that will be used to save the values of the velocity
        field.
    writePath : str
        The path for the writer.
    dt : float
        The time-step to be used in the simulation. This will be used to
        associate a time-value with the produced velocity
        fields.
    t0 : float
        The starting time to be used in the simulation. This will be
        used to associate a time-value with the produced velocity.
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
        The values of eta for the corresponding values of the mean
        velocity from the precursor.
    yPlusPrec : ndarray
        The values of y+ for the corresponding values of the mean
        velocity from the precursor.
    pointsZ : ndarray
        A 2d array containing the values of z for the points of the
        precuror mesh.
    etaInfl : 1d ndarray
        The values of eta for the mesh points at the inflow boundary.
    yPlusInfl : 1d ndarray
        The values of y+ for the meshpoints at the inflow boundary.
    pointsZInfl : int
        A 2d array containing the values of z for the points of the
        inflow boundary.
    nInfl : int
        The amount of points in the wall-normal direction that contain
        the boundary layer at the inflow boundary. That is, for points
        beyound nInfl, 0 will be prescribed.
    nInner : int
        The amount of points where inner rescaling should be considered.
        For points beyond nInner, the outer rescaling only we be
        computed. The relaxes the demand on Re_tau for the precursor.
    gamma : float
        The ration of the friction velocities in the inflow boundary
        layer and the precursor.
    yInd : ndarray
        The sort indices for sorting the read velocity field. This is
        needed when some sorting is performed when the mesh points are
        read and turned into ordered 2d arrays. Them the exact same
        sorting should be applyed to the velocity fields.
    zInd : ndarray
        Same as yInd, but for the sorting of the z values.
    surfaceName : str, optional
        For the foamFile reader, the name of the surface used for
        sampling the velocity values.
    times : list of floats, optional
        For the foamFile reader, the times for which the velocity field
        was sampled in the precursor simulation.
    """

    # Grab info regarding parallelization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProcs = comm.Get_size()

    # Get the total amount of rescalings to be done
    size = int((tEnd-t0)/dt+1)
    
    # Calculate the amount of rescalings each processor is responsible for
    [chunks, offsets] = chunks_and_offsets(nProcs, size)

    # Perform the rescaling
    for i in xrange(chunks[rank]):
        t = t0 + dt*i + dt*int(offsets[rank])
        t = float(("{0:."+str(timePrecision)+"f}").format(t))
        position = int(offsets[rank]) + i

        if rank == 0:
            print "     Rescaled about", i/float(chunks[rank])*100, "%" 

        # Read U data
        if reader == "foamFile":
            assert position < len(times)

            [uPrimeX, uPrimeY, uPrimeZ] = read_u_from_foamfile(
                os.path.join(readPath, times[position], surfaceName,
                             "vectorField", "U"), pointsZ.shape[0],
                pointsZ.shape[1], yInd, zInd, interpolate=True)
        else:
            raise ValueError("Unknown reader")

        # Subtract mean
        uPrimeX -= uMeanPrec[:, np.newaxis]

        [uXInfl, uYInfl, uZInfl] = \
            lund_rescale_fluctuations(
                etaPrec, yPlusPrec, pointsZ,
                uPrimeX, uPrimeY, uPrimeZ, gamma,
                etaInfl, yPlusInfl, pointsZInfl,
                nInfl, nInner)

        # Add mean
        uXInfl += uMeanInfl

        # Write
        if writer == "tvmfv":
            write_velocity_to_tvmfv(writePath, t, uXInfl, uYInfl, uZInfl)
        elif writer == "hdf5":
            write_velocity_to_hdf5(writePath, t, uXInfl, uYInfl, uZInfl,
                                   position)
        else:
            raise ValueError("Unknown writer")
