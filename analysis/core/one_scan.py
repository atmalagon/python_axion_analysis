import numpy as np

from core_utils.features import *
from core_utils.vis_utils import *
from core_utils.fit_utils import *


class Scan(object):
    def __init__(self, param_file, spectrum_file, n_avg=1, C=0.5, beta=1, deg=3):
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
        self.id = params['digitizer_log_reference']

        self.num_averages = n_avg
        self.C = C
        self.beta = beta

        self.freq = np.linspace(self.start, self.stop, self.num_points) #MHz
        self.data = np.array(spectrum['power_spectrum_channel_one']) #Watts

        self.fit = np.array([])
        self.bad_data = 0

    def compute_fit(self, model=simplex_fit, deg=5):
        """
        Choose fit to use from fit_utils: simplex_fit, leastsq_fit,
        poly_fit.
        All models try to fit a polynomial of degree = deg to the data,
        but they have different minimization strategies.
        """

        #(TODO) allow self.uncertainty to be passed in too
        fit = model(self.freq, self.data, deg)
        
        return fit

    def quality_cuts(self):
        """Flags data with poor experimental values."""
        # (TODO) test for isnumeric for params too

        if (scan.Q < 10000) or (scan.B < 1.3) or (scan.Tn > 1.3):
            scan.bad_data = 1

        for arr in [self.freq, self.data]:
            if len(self.freq) != self.num_points:
                scan.bad_data = 1
        
        if np.mean(scan.data) < 1.e-9:
            scan.bad_data = 1

    def get_residuals(self, fit):
        """                                                                               
        Rescale data so that mean is at k_B * Tn * bin_width [Watts].                     
        data is an array of power values in Watts, Tn is a noise temperature              
        in Kelvin, and bin_width is the BW of each bin in Hz.                             
        Return residuals (data - kTB) and uncertainties.
        """
        k_B = 1.3806488e-23 # Watts/Hz/Kelvin                                             

        kbt = (k_B * self.bin_width * 1.e6 * self.Tn)
        sigma = kbt / np.sqrt(self.num_averages)

        #scale data to have mean of 1 in each bin.
        #the fit is approximating the mean in each bin.
        normalize = np.divide(self.data, np.array(fit))

        #residuals are the fluctuations above mean of kbt
        residuals = (normalize - 1) * kbt
        uncertainties = np.full(len(residuals), sigma)

        return residuals, uncertainties

    def create_axion_array(self):
        """
        Creates an array with the axion power * Lorentzian for the same frequencies
        as the scan.
        """
        i = 0
        axion = np.ones(self.num_points)
        for f in self.freq:
            response = lorentz(f, self.mode_f, self.Q)
            on_resonance_pwr =axion_power(self.B, self.mode_f, f - self.mode_f, self.Q)
            axion[i] = response * on_resonance_pwr
            i+=1

        return axion

    def process_and_plot_scan(self, model=simplex_fit, deg=5):
        """
        Cut irrelevant data, scale the power to be at kBT.
        Get fluctuations (power - kBT), and assign uncertainty:
        kBT/sqrt(avgs).
        Finally, scale the power/uncertainties once more to be in units of axion
        power (KSVZ power * Lorentzian)
        """
        self.quality_cuts()

        if self.bad_data == 1:
            print 'data is bad.'
            return

        baseline = self.compute_fit(model, deg)
        residuals, uncertainties = self.get_residuals(baseline)
        axion = self.create_axion_array()

        #plot data with fit
        plot_errorbars(self.freq, self.data, fit=baseline, caption='single_scan_with_fit')

        #plot residuals with their distribution in a hist on the side
        plot_scan_with_hist(self.freq, residuals, label=self.id, caption='single_scan')

        #also show the predicted axion power
        plot_scan_with_hist(self.freq, residuals, prediction=axion, label=self.id,
                            caption='prediction_single_scan')


if __name__ == "__main__":
    path = '../../data/samples/five_scans/'
    param_file = path + 'parameters/one_scan_parameters.npy'
    spectrum_file = path + 'spectra/one_scan.npy'


    scan = Scan(param_file, spectrum_file)

    #hack since all the sample data has broken sensor for squid temperature
    scan.Tn = 0.9

    scan.process_and_plot_scan(model=leastsq_fit, deg=1)
    


