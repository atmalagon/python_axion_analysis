import numpy as np
from scipy import special


def confidence_level(n_sigma, mu, power_sigma):
    """
    Integrate a gaussian from 0 to 1 with mean at n_sigma, to find
    what the confidence level of a particular power level (in units of sigma)
    is.
    """
    summand = 0
    step = 0.002
    for x in np.arange(0., 1., step):
        value = power * np.sinc(x / 2.) ** 2 - n_sigma
        summand += 1 / 2. * (1 + special.erf(value) ) * step
    return summand

def flag_candidates(threshold_sigma, freq, data, error):
    """Flags all data bins with power threshold_sigma * sigma above 0."""
    assert len(data) == len(error)
    candidates_freq = []
    for f, power, fluctuation_size in zip(freq, data, error):
        if power > threshold_sigma * fluctuation_size:
            candidates_freq.append(f)
    return candidates_freq

    
