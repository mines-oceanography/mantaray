import numpy as np

def c_shallow(depth):
    '''
    Shallow water approximation for phase velocity c
    Args:
        depth (float): Depth (positive, measured from sea surface to bathymetry)
    Returns:
        c (float): Phase speed of wave
    '''
    return np.sqrt(9.8*depth)

def c_deep(k):
    '''
    Deep water approximation for phase velocity c
    Args:
        k (float): Wavenumber
    Returns:
        c (float): Phase speed of wave
    '''
    return np.sqrt(9.8/k)

def c_and_cg_calc(depth, k):
    '''
    Calculation of c and cg using full dispersion relationship
    Args:
        depth (float): Depth (positive, measured from sea surface to bathymetry)
        k (float): Wavenumber
    Returns:
        c (float): Phase speed of wave
        cg (float): Group velocity of wave
    '''
    c = np.sqrt(9.8/k * np.tanh(k*depth))
    n = 0.5 * (1 + 2*k*depth / np.sinh(2*k*depth))
    cg = n*c
    return c, cg