import glob
import os

import numpy as np

from one_scan import Scan
from core_utils.vis_utils import plot_overlay


# (TODO) implement this class
# class Overlapped_Scans(Scan):


if __name__ == "__main__":

    #sort by modification time
    param_files = sorted(glob.glob('../../data/samples/ninetynine_scans/parameters/*'), key=os.path.getmtime)
    spectrum_files = sorted(glob.glob('../../data/samples/ninetynine_scans/spectra/*'), key=os.path.getmtime)


    freq_list = []
    array_list = []
    for param_file, spectrum_file in zip(param_files[3:], spectrum_files[3:]):
        scan = Scan(param_file, spectrum_file)
        if scan.mode_f != 696.91:
            array_list.append(scan.data)
            freq_list.append(scan.freq)
            print scan.mode_f

    plot_overlay(freq_list, array_list, start=1, offset=0.7, caption='ninetynine_scans')
    
