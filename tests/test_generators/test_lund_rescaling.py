import eddylicious
from eddylicious.generators.lund_rescaling import *
import numpy as np
import pytest
from os import path
from numpy.testing import assert_almost_equal


#Test recaling dns data onto itself
def test_lund_rescale_mean_velocity():
    file = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                     "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(file)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[: , 1]
    assert_almost_equal(lund_rescale_mean_velocity(eta, yPlus, u, eta.size,
                                                   eta.size, eta, yPlus, 1,
                                                   u[-1], u[-1], 1)[:, 0],
                        u)
