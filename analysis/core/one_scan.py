import numpy as np

from core_utils.features import *
from core_utils.vis_utils import *


class Scan(object):
    def __init__(self, param_file, spectrum_file, n_avg=1, C=0.5, beta=1):
        params = np.load(param_file).item()            
        spectrum = np.load(spectrum_file).item()

        #span is in MHz
        self.start = spectrum['start_frequency_channel_one'] 
        self.stop = spectrum['stop_frequency_channel_one']
        self.span = self.stop - self.start
        self.mode_f = params['mode_frequency_channel_one']

        self.num_points = spectrum['mixing_window_length']
        self.bin_width = self.span / self.num_points # in MHz

        self.Tn = params['cavity_top_temperature'] + params['squid_temperature'] #Kelvin
        self.Q = params['q_channel_one']
        self.B = params['b_field'] #Tesla
        self.digitizer_id = params['digitizer_log_reference']

        self.num_averages = n_avg
        self.C = C
        self.beta = beta

        self.freq = np.linspace(self.start, self.stop, self.num_points) #MHz
        self.data = np.array(spectrum['power_spectrum_channel_one']) #Watts
        # (TODO) make separate code that figures out which parameters minimize residuals.
        self.mean = savitzky_golay(self.data, 11, 4)

        # uncertainty is equal to power over sqrt avgs if the noise is Gaussian
        # (TODO) I'm not sure this is right - is there some sort of pre-averaging
        # going on in the raw data that I need to account for?
        self.uncertainty = self.data / np.sqrt(self.num_averages)

        self.axion = np.zeros(self.num_points)
        self.snr = np.zeros(self.num_points)

        self.bad_data = 0

    def quality_cuts(self):
        """Flags data with poor experimental values."""
        # (TODO) test for isnumeric for params too

        if scan.Q < 10000 | scan.B < 1.3 | scan.Tn > 1.3:
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

        #scale data to have mean of 1. Don't think I need to worry about 
        #div by zero problems?                                                                                 
        normalize = np.divide(data, self.mean)

        rescaled = normalize * (k_B * self.bin_width * 1.e6 * self.Tn)
        return rescaled

    def weight_by_axion_array(self):
        """
        Creates an array with the axion power * Lorentzian for the same frequencies
        as the scan.
        """
        self.axion = np.array([lorentz(f, self.mode_f, self.bin_width) * 
                               axion_power(self.B, self.mode_f, f - self.mode_f,
                                           self.Q, self.Tn) for f in self.freq])

    def process_scan(self):
        """
        Cut irrelevant data, scale the power/uncertainty arrays to be at kBT.
        Finally, scale the power/uncertainties once more to be in units of axion
        power (KSVZ power * Lorentzian)
        """
        if self.bad_data == 1:
            return

        #(TODO) pass in mask limits into process_scan directly?
        self.freq = mask_data(self.freq)

        # (TODO) really need to make rescale_scan a staticmethod
        # (TODO) why doesn't it update when I assign rescaled values to self.data an
        # self.uncertainty?
        rescaled_data = self.rescale_scan(mask_data(self.data))
        rescaled_uncertainty = self.rescale_scan(mask_data(self.uncertainty))

        return rescaled_data, rescaled_uncertainty

if __name__ == "__main__":
    path = '../../data/samples/five_scans/'
    param_file = path + 'parameters/one_scan_parameters.npy'
    spectrum_file = path + 'spectra/one_scan.npy'
    scan = Scan(param_file, spectrum_file)

    # (TDOO) add autoscaling to plotting functions
    plot_errorbars(scan.freq, scan.data, fit=scan.mean, caption='single_scan_with_fit')
    plot_errorbars(scan.freq, scan.data, scan.uncertainty, fit=scan.mean, caption='single_scan')

    rs_data, rs_uncertainty = scan.process_scan()
    plot_errorbars(scan.freq, rs_data, rs_uncertainty, caption='kBT_rescaled_single_scan')
    
    plot_scan_with_hist(scan.freq, rs_data, start=4, caption='single_scan')
