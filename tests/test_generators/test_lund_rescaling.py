import eddylicious
from eddylicious.generators.lund_rescaling import *
from eddylicious.generators.helper_functions import *
import numpy as np
from os import path
from numpy.testing import assert_almost_equal
from scipy.interpolate import interp1d


# Test rescaling dns data onto itself, using the same grid
def test_lund_rescale_mean_velocity_same_grid():
    fileName = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                         "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(fileName)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    assert_almost_equal(lund_rescale_mean_velocity(eta, yPlus, u, eta.size,
                                                   eta.size, eta, yPlus, 1,
                                                   u[-1], u[-1], 1,
                                                   blending_function)[:, 0],
                        u)


# Test rescaling dns data onto itself, using a different grid
def test_lund_rescale_mean_velocity_different_grid():
    fileName = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                         "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(fileName)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    etaInfl = np.linspace(0, eta[-1], 119)
    yPlusInfl = etaInfl*6.37309e-2/3.5e-4
    uTest = interp1d(eta, u)(etaInfl)
    uRescaled = lund_rescale_mean_velocity(eta, yPlus, u, etaInfl.size,
                                           etaInfl.size, etaInfl, yPlusInfl, 1,
                                           u[-1], u[-1], 1,
                                           blending_function)[:, 0]
    for i in range(len(uTest)):
        assert_almost_equal(uTest[i], uRescaled[i], decimal=3)


# Test recaling dns data onto itself, using a different grid with eta > 1
def test_lund_rescale_mean_velocity_eta_greater_one():
    file = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                     "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(file)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    etaInfl = np.linspace(0, 2*eta[-1], 119)
    yPlusInfl = etaInfl*6.37309e-2/3.5e-4

    nInfl = 0
    for i in range(etaInfl.size):
        if etaInfl[i] <= eta[-1]:
            nInfl += 1

    uTest = interp1d(eta, u)(etaInfl[:nInfl])

    uRescaled = lund_rescale_mean_velocity(eta, yPlus, u, nInfl,
                                           nInfl, etaInfl, yPlusInfl, 1,
                                           u[-1], u[-1], 1,
                                           blending_function)[:, 0]
    for i in range(nInfl):
        assert_almost_equal(uTest[i], uRescaled[i], decimal=3)

# Test recaling dns data onto itself, using a different grid with eta > 1
# and nInner < nInfl
def test_lund_rescale_mean_velocity_eta_greater_one_ninner_less_ninfl():
    file = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                     "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(file)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    etaInfl = np.linspace(0, 2*eta[-1], 119)
    yPlusInfl = etaInfl*6.37309e-2/3.5e-4

    nInfl = 0
    for i in range(etaInfl.size):
        if etaInfl[i] <= eta[-1]:
            nInfl += 1

    uTest = interp1d(eta, u)(etaInfl[:nInfl])

    uRescaled = lund_rescale_mean_velocity(eta, yPlus, u, nInfl,
                                           int(0.5*nInfl), etaInfl, yPlusInfl, 1,
                                           u[-1], u[-1], 1,
                                           blending_function)[:, 0]
    for i in range(nInfl):
        assert_almost_equal(uTest[i], uRescaled[i], decimal=3)
