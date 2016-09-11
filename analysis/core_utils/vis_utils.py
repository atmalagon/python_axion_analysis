import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import norm


matplotlib.style.use('ggplot')
plt.rcParams['figure.figsize'] = (8, 6)
ticklabel_format(style='sci', axis='x', scilimits=(0,0))

pltdir = '../../plots/'

def plot_errorbars(X, Y, error, LO, caption='', additional_offset=4.09):
    """
    Saves 2 plots to file:
    (1) Y vs X with error bars for positive X only,
    (2) Y vs X with error bars for all X.
    The frequency axis is shown relative to the LO frequency.
    The additional_offset is a relic from the YMCE analysis.
    """
    # (TODO) saving the positive frequencies only is a relic of the YMCE analysis.
    # Remove for ADMX analysis?
    id = np.where(X>0.)[0][0]
    for start in [0, id]:
        plt.errorbar(X[start:]/1.e6, Y[start:], yerr=error[start:], fmt='o')
        plt.xlabel('Freq Offset (MHz) from {0} GHz'.format(str(LO + additional_offset)))
        plt.ylabel('$\delta P$ [mW/Hz]')
        plt.grid(True)
        plt.title(caption)
        plt.savefig(pltdir + caption + '_start_' + start + '.png', bbox_inches='tight')
        plt.close()


def plot_hist(pow_hist, set_bins=30, caption=''):
    """Plots and saves  a histogram of the input list 'hist'. """
    (mu, sigma) = norm.fit(pow_hist)
    n, bins, patches = hist(pow_hist, bins=set_bins, color='b', normed=1, histtype='step',
                            label='data')
    hist_fit = normpdf(bins, mu, sigma)
    plt.plot(bins, hist_fit, 'r', label='fit')
    plt.xlabel('$\delta P$ [mW/Hz]')
    plt.ylabel('Count')
    plt.legend()
    plt.title(caption)
    plt.savefig(pltdir + caption + '_histogram.png', bbox_inches='tight')
    plt.close()

def plot_coadded(coadded_freq, coadded_pow, coadded_weight, err, caption=''):
    """Saves plot of coadded weighted spectrum. See docs for details on the methodology."""
    weighted_pow = [np.average(p, weights=w) for p, w in zip(coadded_pow, coadded_weight)]
    
    plt.errorbar(coadded_freq, weighted_pow, yerr=
