import glob
import os

def load_test_files():
    path = '~/python_axion_analysis/data/samples/ninetynine_scans/'
    param_files = sorted(glob.glob(path+'/parameters/*'), key=os.path.getmtime)
    spectrum_files = sorted(glob.glob(path+'/spectra/*'), key=os.path.getmtime)

    return param_files, spectrum_files
