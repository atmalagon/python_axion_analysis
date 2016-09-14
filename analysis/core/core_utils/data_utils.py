import glob
import os

def load_test_files(path):
    """
    Loads the list of parameter and spectrum files in the path..
    """
    param_files = sorted(glob.glob(path+'/parameters/*.npy'), key=os.path.getmtime)
    spectrum_files = sorted(glob.glob(path+'/spectra/*.npy'), key=os.path.getmtime)
    return param_files, spectrum_files

if __name__=="__main__":
    p, s = load_test_files('../../../data/samples/ninetynine_scans')
    assert len(p) > 0
    assert len(p) == len(s)
