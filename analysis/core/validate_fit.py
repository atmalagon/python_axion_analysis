import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import interp1d
from scipy.optimize import brute

from core_utils.features import savitzky_golay, baseline_als
from core_utils.data_utils import load_test_files
from one_scan import Scan

def plot_als_tests():
    """
    Trying to use a gridsearch function to optimize
    paramaters for asymmetric least squares smoothing function.
    """
    path = '../../data/samples/ninetynine_scans'
    param_files, spectrum_files = load_test_files(path)
    scan = Scan(param_files[3], spectrum_files[3])

    lam_list = slice(50, 2400, 100) 
    p_list = slice(0.1, 0.9, .05)
    def obj_func(z):
        lam, p = z
        return scan.data - baseline_als(scan.data, lam, p)

    x0, fval,grid ,Jout = brute(obj_func, (lam_list, p_list),full_output=True)

    print x0
    print fval
    print grid
    print Jout
    return

def plot_savgol_tests():
    """
    My attempt to make a training dataset to see how the rms error
    of the savitzky_golay smooth changed with a varying window size.
    Followed:
    http://www.astroml.org/book_figures/chapter8/
    fig_cross_val_C.html
    """
    param_files, spectrum_files = load_data()
    windows = np.arange(7, 55, 2)

    fig = plt.figure(figsize=(8, 6))

    for order, subplot in zip([2, 3, 4], [311, 312, 313]):
        ax = fig.add_subplot(subplot)
        train_err = np.zeros(windows.shape)

        scan = Scan(param_files[3], spectrum_files[3])        
        #not sure I understood this, but the point was
        #to take new data points and compare them to what
        #the filter thought should be the case, based on
        #the data it'd gotten.

        y_interp = interp1d(scan.freq, scan.data)
        x_lin = np.linspace(np.min(scan.freq), np.max(scan.freq), 10000)
        mask = np.random.choice(x_lin, 250, replace=False)
        y_train = y_interp(mask)

        for i, window in enumerate(windows):
            filter = savitzky_golay(y, window, order)
            train_err[i] = np.sqrt(np.sum( (filter[window:] - y[window:]) ** 2) 
                  / len(y[window:]))

        #ax.plot(scan.freq[window:], filter[window:])
        #ax.plot(scan.freq[window:], y_train[window:])
        ax.plot(windows, train_err, label='order {}'.format(order))

        ax.legend(loc=2)
        plt.setp(ax.get_xticklabels()[:1], visible=False)
        plt.setp(ax.get_xticklines()[:1], visible=False)

    plt.xlabel('window size of filter')
    plt.ylabel('rms error')
    plt.show()

if __name__=="__main__":
    plot_als_tests()
