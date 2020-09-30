# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

"""Module containing functions for the generation of inflow fields.

"""
from .lund_rescaling import *
from .interpolation import *
from .sem import *
from .random_fluctuations import *

__all__ = ["lund_rescaling", "interpolation", "sem", "random_fluctuations"]
__all__.extend(lund_rescaling.__all__)
__all__.extend(interpolation.__all__)
__all__.extend(sem.__all__)
__all__.extend(random_fluctuations.__all__)

