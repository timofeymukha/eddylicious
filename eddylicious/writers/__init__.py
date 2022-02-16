# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Module containing functions for writing out the inflow fields in various
file formats.

"""

from .hdf5_writers import *
from .ofnative_writers import *
from .vtk_writers import *
from .ram_writers import *

__all__ = ["hdf5_writers", "ofnative_writers", "vtk_writers", "ram_writers"]
__all__.extend(hdf5_writers.__all__)
__all__.extend(ofnative_writers.__all__)
__all__.extend(vtk_writers.__all__)
__all__.extend(ram_writers.__all__)
