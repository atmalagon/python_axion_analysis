import numpy as np
import scipy as sp

from math import factorial

def baseline_als(y, lam, p, niter=10):
    """
    Asymmetric least squares smoothing. From
    https://stackoverflow.com/questions/29156532/
    python-baseline-correction-library.
    y is the data, lam (for lambda) is the weight given
    for how strong you want smoothing to be (usually 100
    to 1.e9); p says how strongly you care about asymmetry 
    (generally 0.001 to 0.1). niter is the number of iterations,
    computing new weights.
    """
    L = len(y)
    D = sp.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in xrange(niter):
        W = sp.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sp.sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
    From the Scipy cookbook: https://scipy.github.io/old-wiki/
    pages/Cookbook/SavitzkyGolay. There is now a scipy.signal.savgol_filter
    implementation. This is a way to empirically get the baseline of the data
    as opposed to a polynomial fit, which sometimes fails for reasons I
    don't understand.
    """
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

def lorentz(f, f_0, Q):
    """Returns a lorentzian centered at f_0 with quality Q."""
    FWHM = f_0 / float(Q)
    return 1. / (1. + 4. * (f - f_0)**2 / FWHM**2)

def mask_data(X, start=0, left=0, right=0, end=256):
    """Returns two columns of data cut on the wings and in the middle."""
    Xmask = np.concatenate((X[start:left], X[right:end]))
    return Xmask

def g_a2gamma_ksvz(f):
    """
    Tells you what g_a2gamma is at KSVZ, in GeV^-1 for
    a given frequency (input as Hz).
    """
    g = 0.97
    alpha = 1. / 137.
    m_times_f = 0.006 #GeV^2
    hbar = 6.582122e-25 # GeV-sec                                                  c = 3.e10 # cm / sec                                                       
    m = hbar * 2 * np.pi * f
    
    g_a2gamma = g * alpha / (np.pi * m_times_f) * m
    return g_a2gamma

def axion_power(B, f_0, delta, Q, C=0.5, beta=1., g=0.97, rho=0.45,
                V=0.0015):
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
    c = 3.e10 # cm / sec

#    V = 0.0015 # volume in m**3 (for ADMX 2014-16)
    omega = 2*np.pi*f_0*1.e6 # radians/sec assuming freq is given in MHz
    FWHM = f_0 / float(Q)

    phenom_terms = (g * alpha / (np.pi * m_times_f))**2 * rho * (hbar * c)**3
    physical_terms = omega / mu * B**2 * V * C * Q * beta / (1. + beta) 
    off_resonance = 1. / (1. + (2. * delta / FWHM)**2)
    axion_power = phenom_terms * physical_terms * off_resonance
    return axion_power # in watts

if __name__ == "__main__":
    #check that axion_power is giving benchmark results
    print axion_power(9, 5000, 0, 12000, beta=1.4)
    #should be 5.e-24

    #check that lorentzian at peak is 1
    print lorentz(5000, 5000, 12000)

    #check that g_a2gamma at 1GHz is ~ 1.5e-15 GeV^-1.
    print g_a2gamma_ksvz(0.7e9)
