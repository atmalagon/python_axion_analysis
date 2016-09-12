import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

from scipy.stats import norm


matplotlib.style.use('ggplot')
plt.rcParams['figure.figsize'] = (8, 6)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

pltdir = '../../plots/'

def plot_errorbars(X, Y, error=None, fit=None, start=0, caption=''):
    """
    Saves plot to disk:
    (1) Y vs X with error bars (optional) and fit superimposed (optional).
    Can also choose starting point of data to plot (starts from beginning of
    input array by default.)
    """
    if error is not None:
        plt.errorbar(X[start:], Y[start:], yerr=error[start:], fmt='o')
    else:
        plt.scatter(X[start:], Y[start:])

    if fit is not None:
        plt.plot(X[start:], fit[start:], label='fit')
        plt.legend()

    plt.xlabel('Frequency [MHz]')
    plt.ylabel('$\delta P$ [Watts]')
    plt.grid(True)
    plt.title(caption.replace('_', ' '))

    plt.savefig(pltdir + caption + '_start_' + str(start) + '.png', bbox_inches='tight')
    plt.close()


def plot_hist(pow_hist, set_bins=30, caption=''):
    """Plots and saves  a histogram of the input list 'hist'. """
    (mu, sigma) = norm.fit(pow_hist)
    n, bins, patches = hist(pow_hist, bins=set_bins, color='b', normed=1, histtype='step',
                            label='data')
    hist_fit = normpdf(bins, mu, sigma)
    plt.plot(bins, hist_fit, 'r', label='fit')
    plt.xlabel('$\delta P$ [Watts]')
    plt.ylabel('Count')
    plt.legend()
    plt.title(caption)
    plt.savefig(pltdir + caption + '_histogram.png', bbox_inches='tight')
    plt.close()

def plot_overlay(freq_list, array_list, start=0, offset=1, caption=''):
    """Plots multiple single scans, with some fixed vertical offset."""
    lift = 0
    min = np.inf
    max = -np.inf
    freq_left = np.inf
    freq_right = -np.inf
    assert len(freq_list) == len(array_list)
    colors = cm.rainbow(np.linspace(0, 1, len(array_list)))
    frame = plt.gca()

    avg_value = np.mean(array_list[0])
    for freq, spectrum, c in zip(freq_list, array_list, colors):
        plt.scatter(freq[start:], spectrum[start:] + lift, color=c)
        lift += offset * avg_value # go up by some percent of the avg pwr
        
    plt.xlim([np.amin(freq_list), np.amax(freq_list)])
    plt.ylim([np.amin(array_list), np.amax(array_list) + lift])
    plt.xlabel('Frequency [MHz]')

    #plt.grid(True)
    frame.axes.get_yaxis().set_visible(False)
    #plt.title(caption.replace('_', ' '))

    plt.savefig(pltdir + 'overlay_' + caption + '_start_' + str(start) + '.png',
                bbox_inches='tight')
    plt.close()
