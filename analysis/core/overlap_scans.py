import glob
import os

import numpy as np

from core_utils.vis_utils import plot_overlay, plot_gsquared
from core_utils.data_utils import load_test_files
from one_scan import Scan


class Overlapped_Scans(object):
    """"
    This class stitches together the gsquared measurements
    from multiple scans, where the stitching is done by
    a weighted average. The weights are 1/err^2. I followed
    Gray's procedure for this.
    """
    def __init__(self, start, stop, bin_width):
        self.start = start
        self.bin_width = bin_width
        self.all_freqs = np.arange(start, stop, bin_width)
        bins = len(self.all_freqs)
        self.total = np.zeros(bins)
        self.weights = np.zeros(bins)

    def walk_through_files(self, param_files, spectrum_files):
        for p, s in zip(param_files, spectrum_files):
            scan = Scan(p, s)
            scan.Tn = 0.66 #hack because test data had wonky sensor reading

            result = scan.process_and_plot_scan(plots_on=False)
            if not result:
                #go to next scan if this scan didn't produce a result.
                continue

            gsquared, err = result
            assert len(gsquared) == len(err)

            index = (scan.start - self.start) // self.bin_width

            for j in range(len(gsquared)):
                idx = int(index) + j
                if self.weights[idx] == 0:
                    #if there was nothing there, fill in
                    #value for gsquared and weight.
                    self.total[idx] = gsquared[j]
                    self.weights[idx] = 1 / err[j]**2
                else:
                    #do a weighted sum, with weights given by
                    #1/err^2, and normalized by sum of weights.
                    temp = self.weights[idx] * self.total[idx]
                    self.weights[idx] += 1 / err[j]**2
                    self.total[idx] = (temp + gsquared[j] / err[j]**2)
                    self.total[idx] /= self.weights[idx]
        
    def plot_total(self, caption):
        """Plot the output gsquared values.."""
        idx = np.arange(20, len(self.all_freqs) - 10)
        plot_gsquared(self.all_freqs[idx], self.total[idx], caption=caption)
        
        candidates = []
        position = []
        mask = []
        for i, val in enumerate(self.total):
            if val > 5 * self.weights[i]:
                candidates.append(val)
                position.append(i)
                mask.append(False)
            else:
                mask.append(True)

        mask = np.array(mask)
        print len(outliers)
        plot_gsquared(self.all_freqs[mask], self.total[mask], caption='total_candidates_removed')
        return position, candidates

if __name__ == "__main__":
    path = '../../data/eng_run'
    p, s = load_test_files(path)
    scan_first = Scan(p[0], s[0])

    combined = Overlapped_Scans(696.69278879272, 697.71158285522, scan_first.bin_width)
    combined.walk_through_files(p, s)
    position, candidates = combined.plot_total(caption='eng_run_gsquared')
    with open('../../data/candidates.txt', 'w') as file:
        print >> file, position, candidates
    print len(position)
#    path = '../../data/samples/ninetynine_scans'
#    freq_list = []
#    array_list = []
#    for param_file, spectrum_file in zip(param_files[3:], spectrum_files[3:]):
#        scan = Scan(param_file, spectrum_file)
#        if scan.mode_f != 696.91:
#            array_list.append(scan.data)
#            freq_list.append(scan.freq)

#    plot_overlay(freq_list, array_list, start=1, offset=0.7, caption='ninetynine_scans')
    
