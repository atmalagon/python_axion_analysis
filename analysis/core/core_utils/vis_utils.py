import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import norm


matplotlib.style.use('ggplot')
plt.rcParams['figure.figsize'] = (8, 6)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

pltdir = '../../plots/'

def plot_errorbars(X, Y, error, fit=None, start=0, caption=''):
    """
    Saves plot to disk:
    (1) Y vs X with error bars for all X.
    """
    plt.errorbar(X[start:], Y[start:], yerr=error[start:], fmt='o')
    if fit is not None:
        plt.plot(X[start:], fit)

    plt.xlabel('Frequency [MHz]')
    plt.ylabel('$\delta P$ [Watts]')
    plt.grid(True)
    plt.title(caption)
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

def plot_coadded(coadded_freq, coadded_pow, coadded_weight, err, caption=''):
    """
    Saves plot of coadded weighted spectrum.
    See docs for details on the methodology.
    """
    weighted_pow = [np.average(p, weights=w) for p, w in zip(coadded_pow,
                                                             coadded_weight)]
