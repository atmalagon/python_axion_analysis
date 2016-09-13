import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import interp1d

from core_utils.features import savitzky_golay
from one_scan import Scan


#let's optimize the savitzky_golay filter on the data
def calculate_filter_rms(y, window, order=3):
    """
    Calculate the rms error of the filter from the data
    for different window sizes, using randomly selected 
    points in the spectrum as training data. Follows:
    http://www.astroml.org/book_figures/chapter8/
    fig_cross_val_C.html
    """
    filter = savitzky_golay(y, window, order)
    err = np.sqrt(np.sum( (filter[window:] - y[window:]) ** 2) 
                  / len(y[window:]))

    return filter, err

path = '../../data/samples/ninetynine_scans/'
param_files = sorted(glob.glob(path+'/parameters/*'), key=os.path.getmtime)
spectrum_files = sorted(glob.glob(path+'/spectra/*'), key=os.path.getmtime)

windows = np.arange(7, 55, 2)

fig = plt.figure(figsize=(8, 6))

for order, subplot in zip([2, 3, 4], [311, 312, 313]):
    ax = fig.add_subplot(subplot)
    train_err = np.zeros(windows.shape)

    scan = Scan(param_files[3], spectrum_files[3])        
    y_interp = interp1d(scan.freq, scan.data)
    x_lin = np.linspace(np.min(scan.freq), np.max(scan.freq), 10000)
    mask = np.random.choice(x_lin, 250, replace=False)
    y_train = y_interp(mask)

    for i, window in enumerate(windows):
        filter, train_err[i] = calculate_filter_rms(y_train, window, order)

        intrinsic_err = [np.std(scan.data[:10])/np.sqrt(10)]*len(windows)
    import ipdb
    ipdb.set_trace()
    ax.plot(scan.freq[window:], filter[window:])
    ax.plot(scan.freq[window:], y_train[window:])
    #ax.plot(windows, train_err, label='order {}'.format(order))
    #ax.plot(windows, intrinsic_err)
    #ax.legend(loc=2)
    plt.setp(ax.get_xticklabels()[:1], visible=False)
    plt.setp(ax.get_xticklines()[:1], visible=False)

#plt.xlabel('window size of filter')
#plt.ylabel('rms error')


plt.show()
