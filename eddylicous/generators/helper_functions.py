import numpy as np
from scipy.integrate import simps


# Blending function for inner and outer scales
def blending_function(eta, alpha=4, b=0.2):
    """Return the value of the blending function W for Lund's rescaling.

    Return the value of the blending function for inner and
    and outer profiles produced by Lund's rescaling.
    For eta>1 the function returns 1.

    Keyword arguments:
    alpha -- the value of alpha (delfault 4.0)
    b -- the value of b (defaul 0.2)
    """

    val = np.zeros(eta.shape)

    for i in xrange(eta.size):
        if eta[i] <= 1:
            val[i] = 0.5*(1+1/np.tanh(alpha)*np.tanh(alpha*(eta[i]-b) /
                          ((1-2*b)*eta[i]+b)))
        else:
            val[i] = 1.0
    return val


def delta_99(y, v):
    """Compute delta_99."""
    for i in xrange(y.size):
        if v[i] >= 0.99*np.max(v):
            delta99 = y[i]
            break
    return delta99


def delta_star(y, v):
    """Compute delta^* using Simpson's method."""
    return simps((1-v/v[-1]), x=y)


def theta(y, v):
    """Compute the momentum thickness using Simpson's method."""
    return simps(v/v[-1]*(1-v/v[-1]), x=y)
