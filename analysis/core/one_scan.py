import numpy as np

from core_utils.feature import *

class Scan(object):
    def __init__(self, param_file, spectrum_file, n_avg=1, C=0.5, beta=1):
        params = np.load(param_file).item()            
        spectrum = np.load(spectrum_file).item()

        #span is in MHz
        self.start = spectrum['start_frequency_channel_one'] 
        stop = spectrum['stop_frequency_channel_one']
        self.span = stop - start

        self.num_points = spectrum['mixing_window_length']
        self.bin_width = span * 1.e6 / scan.num_points # in Hz

        self.Tn = params['cavity_top_temperature'] + params['squid_temperature'] #Kelvin
        self.Q = params['q_channel_one']
        self.B = params['b_field'] #Tesla
        self.digitizer_id = params['digitizer_log_reference']

        self.num_averages = n_avg
        self.C = C
        self.beta=beta

        self.freq = np.arange(scan.start, scan.start + scan.span, scan.bin_width)
        self.data = np.array(spectrum['power_spectrum_channel_one'])
        #uncertainty is equal to power over sqrt avgs if the noise is Gaussian
        self.uncertainty = scan.data / np.sqrt(scan.num_averages)

        self.bad_data = 0

    def quality_cuts(self):
        """Flags data with poor experimental values."""
        # (TODO) test for isnumeric for params too

        if scan.Q < 10000 | scan.B < 1.3 | scan.Tn > 1:
            scan.bad_data = 1

        for arr in [self.freq, self.data]:
            if len(self.freq) != self.num_points:
                scan.bad_data = 1

    def rescale_scan(self, data):
        """                                                                               
        Rescale data so that mean is at k_B * Tn * bin_width [Watts].                     
        data is an array of power values in Watts, Tn is a noise temperature              
        in Kelvin, and bin_width is the BW of each bin in Hz.                             
        """
        #(TODO) should I be using a staticmethod here?

        k_B = 1.3806488e-23 # Watts/Hz/Kelvin                                             

        #(TODO) make separate code that figures out which parameters minimize residuals.
        mean = savitzky_golay(data, 11, 4)
        #(TODO): pass in smoothing parameters as input to main function?                  

        #scale data to have mean of 1. Don't think I need to worry about 
        #div by zero problems?                                                                                 
        normalize = np.divide(data, mean)

        rescaled = np.multiply(normalize, k_B * self.bin_width * self.Tn)
        return rescaled

    def process_scan(self):
        """
        Cut irrelevant data, scale the power/uncertainty arrays.
        Finally, scale the power/uncertainties once more to be in units of axion
        power (KSVZ power * Lorentzian)
        """
        if self.bad_data == 1:
            return

        #(TODO) pass in mask limits into process_scan directly?
        for arr in [self.freq, self.data, self.uncertainty]:
            arr = rescale_scan(mask_data(arr))
            

