# This file is part of eddylicious
# (c) 2020 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for generating inflow velocity fields using
the Synthetic Eddy Method

"""
from __future__ import print_function
from __future__ import division
import numpy as np
from mpi4py import MPI
from ..helper_functions import chunks_and_offsets
from ..writers.ofnative_writers import write_velocity_to_ofnative
from ..writers.hdf5_writers import write_velocity_to_hdf5

__all__ = ["sem_generate"]

def sem_generate(inflowPoints,
                 sigma,
                 reynoldsStress):
    """Generate the the inflow velocity using Lund's
    rescaling.

    This function will use Lund et al's rescaling in order to generate

    Parameters
    ----------
    readerFunction : function
        The function to use for reading in data, generated by the
        reader. Should contain the reader's name in the attribute
        "reader".
    """

    # Determine the properties of the box to be filled with eddies

    boxMin = np.min(inflowPoints - sigma, 0)
    boxMax = np.max(inflowPoints - sigma, 0)
    boxVolume = np.prod(boxMax - boxMin)




    pass