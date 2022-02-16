# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Functions for reading fields stored in generic binary files.

"""
import numpy as np
import os
from ..helper_functions import read_float_array

__all__ = ["read_structured_velocity_binary"]


def read_structured_velocity_binary(readPath: str, eMode: str, wordSize: int,
                                    ny: int, nz: int, order: str):
    """Read in an float array of size ny x nz from a binary file.

    Parameters
    ----------
    readPath : str
        Path to the file
    eMode : {'<', '>'}
        Endian mode.
    wordSize : {4, 8}
        Size of the float, corresponding to single or double precision.
    ny : int
        Number of values in the y direction.
    nz : int
        Number of values in the z direction.
    order : {'C', 'F'}
        C or Fortran order when reshaping the data to ny x nz dimenstions.

    Returns
    -------
    ndarray
        The read in values

    """
    file = open(readPath, 'rb')
    return read_float_array(file, eMode, wordSize, ny * nz).reshape(ny, nz,
                                                                    order=order)


