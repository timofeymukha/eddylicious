# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

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
    v = np.zeros(u.shape)
    w = blending_function(eta)
    assert_almost_equal(lund_rescale_mean_velocity(eta, yPlus, u, v, eta.size,
                                                   eta, yPlus, 1,
                                                   u[-1], u[-1], 1,
                                                   w)[0][:,0],
                       u)


# Test rescaling dns data onto itself, using a different grid
def test_lund_rescale_mean_velocity_different_grid():
    fileName = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                         "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(fileName)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    v = np.zeros(u.shape)
    etaInfl = np.linspace(0, eta[-1], 119)
    yPlusInfl = etaInfl*6.37309e-2/3.5e-4
    uTest = interp1d(eta, u)(etaInfl)
    w = blending_function(etaInfl)
    uRescaled = lund_rescale_mean_velocity(eta, yPlus, u, v, etaInfl.size,
                                           etaInfl, yPlusInfl, 1,
                                           u[-1], u[-1], 1,
                                           w)[0][:, 0]
    for i in range(len(uTest)):
        assert_almost_equal(uTest[i], uRescaled[i], decimal=3)


# Test recaling dns data onto itself, using a different grid with eta > 1
def test_lund_rescale_mean_velocity_eta_greater_one():
    fileName = path.join(eddylicious.__path__[0], "..", "tests", "datasets",
                         "channel_flow_180", "dns.dat")
    dns = np.genfromtxt(fileName)

    eta = dns[:, 0]
    yPlus = dns[:, 1]
    u = dns[:, 1]
    etaInfl = np.linspace(0, 2*eta[-1], 119)
    yPlusInfl = etaInfl*6.37309e-2/3.5e-4
    v = np.zeros(u.shape)
    w = blending_function(etaInfl)

    nInfl = 0
    for i in range(etaInfl.size):
        if etaInfl[i] <= eta[-1]:
            nInfl += 1

    uTest = interp1d(eta, u)(etaInfl[:nInfl])

    uRescaled = lund_rescale_mean_velocity(eta, yPlus, u, v, nInfl,
                                           etaInfl, yPlusInfl, 1,
                                           u[-1], u[-1], 1,
                                           w)[0][:, 0]
    for i in range(nInfl):
        assert_almost_equal(uTest[i], uRescaled[i], decimal=3)

