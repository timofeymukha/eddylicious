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
            delta99 = y[i-1]
            break
    return delta99


def delta_star(y, v):
    """Compute delta^* using Simpson's method."""
    return simps((1-v/v[-1]), x=y)


def theta(y, v):
    """Compute the momentum thickness using Simpson's method."""
    return simps(v/v[-1]*(1-v/v[-1]), x=y)

def chunks_and_offsets(nProcs, size):
    """Given the size of a 1d array and the number of
    processors, compute chunk-sizes for each processor
    and the starting indices (offsets) for each 
    processor.

    Parameters
    ----------
        nProcs : int
            The amount of processors.
        size : int
            The size of the 1d array to be distributed.

    Returns
    -------
        List of two 1d ndarrays of size nProcs
        The first array contains the chunk-size for each
        processor.
        The second array contains the offset (starting 
        index) for each processor.
    """

    
    chunks = np.zeros((nProcs, 1), dtype=np.int64)
    nrAlloced = 0
    for i in xrange(nProcs):
        remainder = size - nrAlloced
        buckets = (nProcs - i)
        chunks[i] = remainder / buckets
        nrAlloced += chunks[i]

    # Calculate the offset for each processor
    offsets = np.zeros(chunks.shape, dtype=np.int64)

    for i in xrange(offsets.shape[0]-1):
        offsets[i+1] = np.sum(chunks[:i+1])

    if (np.sum(chunks) != size):
        print "Apacha!"
            
    return [chunks, offsets]

def chauhan_U_inner(yPlus, kappa=0.384, a=-10.361):
    alpha = (-1/kappa - a)/2
    beta = np.sqrt(-2*a*alpha-alpha**2)
    R = np.sqrt(alpha**2 + beta**2)
    return  1/kappa*np.log(-(yPlus-a)/a) + R**2/(a*(4*alpha-a))* \
            ((4*alpha+a)*np.log(-a/R*np.sqrt((yPlus-alpha)**2 + beta**2)/ \
            (yPlus - a)) + alpha/beta*(4*alpha + 5*a)*(np.arctan((yPlus-alpha)/beta) \
            + np.arctan(alpha/beta)))

def chauhan_U_inner_mod(yPlus, kappa=0.384, a=-10.361):
    return chauhan_U_inner(yPlus, kappa, a) + 1/2.85*np.exp(-np.log(yPlus/30)**2)

def chauhan_wake(eta, Pi, a2=132.8410, a3=-166.2041, a4=71.9114):
    nom1 = 1 - np.exp(-0.25*(5*a2 + 6*a3 + 7*a4)*eta**4 + a2*eta**5 + a3*eta**6 + a4*eta**7)
    nom2 = 1 - 0.5/Pi*np.log(eta)
    denom = 1 - np.exp(-0.25*(a2+2*a3+3*a4))
    return nom1*nom2/denom

def chauhan_U_composite(yPlus, eta, Pi, kappa=0.384, a=-10.361):
    return chauhan_U_inner_mod(yPlus, kappa, a) + 2*Pi/kappa*chauhan_wake(eta, Pi)

def epsilon_ReT(y, ReTau, kappa, APlus):
    return (0.5*np.sqrt(1 + kappa**2*ReTau**2/9*(1 - (y - 1)**2)**2*(1 + 2*(y-1)**2)**2*(1 - np.exp(-y*ReTau/APlus))**2) - 0.5)

