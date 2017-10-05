# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Module containing functions for reading in data from various file formats.

"""
from .foamfile_readers import *
from .hdf5_readers import *

__all__ = ["foamfile_readers", "hdf5_readers"]
__all__.extend(foamfile_readers.__all__)
__all__.extend(hdf5_readers.__all__)
