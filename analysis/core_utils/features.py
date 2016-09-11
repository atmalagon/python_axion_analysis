import numpy as np
from math import factorial

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def lorentz(f, f_0, BW):
    """Returns a lorentzian centered at f_0 with FWHM of  BW."""
    return 1. / (1. + 4. * (f - f_0)**2 / BW**2)

def mask_data(X, Y, start=30, left=137, right=157, end=264):
    """Returns two columns of data cut on the wings and in the middle."""
    assert len(X) == len(Y)
    assert len(X) > end

    #(TODO) make the defaults more specific to ADMX

    Xmask = np.concatenate((X[start:left], X[right:end]))
    Ymask = np.concatenate((Y[start:left], Y[right:end]))
    return Xmask, Ymask

def axion_power(B, f, delta, Q, Tnoise, C=0.5, beta=1, g=0.97, rho=0.45):
    """
    Returns the expected axion power for the given experimental parameters.
    Taken following Gp Carosi's code, which in turn followed Ben Brubaker's
    writeup. g is the dimensionless axion coupling strength for different
    models; the default is for the KSVZ model. rho is the estimated dark
    matter density in our Milky Way in GeV/cm**3. C is the cavity form
    factor and beta the antenna coupling (both unitless).
    """
    alpha = 1. / 137. # fine structure constant
    m_times_f = 0.006 # axion mass times scale constant in GeV**2
    hbar = 6.582122e-25 # GeV-sec
    mu = np.pi * 4.e-7 # permeability in kg m/charge**2

    V = 0.0015 # volume in m**3 (for ADMX 2014-16)
    omega = 2*np.pi*f*1.e6 # radians/sec assuming freq is given in MHz

    phenom_terms = (g * alpha / (np.pi * m_times_f))**2 * rho * (hbar * c)**3
    physical_terms = omega / mu * B**2 * V * C * Q * beta / (1 + beta) 
    off_resonance = 1 / (1 + (2 * delta)**2)
    axion_power = phenom_terms * physical_terms * off_resonance
    return axion_power # in watts
