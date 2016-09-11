import numpy as np

from core_utils.feature import *


def process_scan(filename):
    """
    Read in a file with a power spectra. Cut irrelevant data, assign an
    array for the uncertainties, and scale the power/uncertainty arrays.
    Finally, scale the power/uncertainties once more to be in units of axion
    power.
    """
    scan = Scan()
    scan.num_averages = 1
    scan.bin_width = xxx
    scan.Tn = xxx

    scan.freq, scan.data = mask_data(scan.freq, scan.data, start=0, left=0,
                                     right=0, end=256)

    #uncertainty is equal to power over sqrt avgs if the noise is Gaussian
    scan.uncertainty = scan.data / np.sqrt(scan.num_averages)

    # (TODO) make rescale_scan a method of Scan the class.
    scan.data = rescale_scan(scan.data, scan.Tn, scan.bin_width)
    scan.uncertainty = rescale_scan(scan.data, scan.Tn, scan.bin_width)

    # [WIP] leave it for now
