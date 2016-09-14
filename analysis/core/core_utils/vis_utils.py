import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from scipy.stats import norm


matplotlib.style.use('ggplot')
plt.rcParams['figure.figsize'] = (8, 6)
plt.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2, 2))

pltdir = '../../plots/'

def plot_errorbars(X, Y, error=None, fit=None, label=None, caption=''):
    """
    Saves plot to disk:
    (1) Y vs X with error bars (optional) and fit superimposed (optional).
    """
    if np.mean(Y) < 1.e-20:
        Yplot = Y * 1.e24
        plt.ylabel('Power [yoctoWatts]')
    else:
        Yplot = Y
        plt.ylabel('Power [Watts]')

    if error is not None:
        plt.errorbar(X, Yplot, yerr=error,label=label)
    else:
        plt.scatter(X, Yplot, label=label)

    if fit is not None:
        plt.plot(X, fit, label='fit')
    
    plt.legend(scatterpoints=1)
    plt.xlabel('Frequency [MHz]')

    plt.grid(True)
    plt.title(caption.replace('_', ' '))

    plt.ylim([0.97* np.amin(Yplot), 1.03 * np.amax(Yplot)])
    plt.savefig(pltdir + caption + '.png', bbox_inches='tight')
    plt.close()

def plot_gsquared(x, gsquared, label=None, caption=''):
    """
    Saves plot of  gsqured vs x, gsquared is in inverse GeV^2.
    """
    plt.scatter(x, gsquared, label=None)

    plt.xlabel('Frequency [MHz]')
    plt.ylabel('$g_{a\gamma\gamma}^2$ [1/GeV^2]')
    plt.grid(True)
    plt.title(caption.replace('_', ' '))

    plt.ylim([0.97* np.amin(gsquared), 1.03 * np.amax(gsquared)])
    plt.savefig(pltdir + caption + '.png', bbox_inches='tight')
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

def plot_overlay(freq_list, array_list, offset=0.7, caption=''):
    """Plots multiple single scans, with some fixed vertical offset."""
    lift = 0
    min = np.inf
    max = -np.inf
    freq_left = np.inf
    freq_right = -np.inf
    assert len(freq_list) == len(array_list)
    colors = cm.jet(np.linspace(0, 1, len(array_list)))
    frame = plt.gca()

    avg_value = np.mean(array_list[0])
    for freq, spectrum, c in zip(freq_list, array_list, colors):
        plt.scatter(freq, spectrum + lift, color=c)
        lift += offset * avg_value # go up by some percent of the avg pwr
        
    plt.xlim([np.amin(freq_list), np.amax(freq_list)])
    plt.ylim([np.amin(array_list), np.amax(array_list) + lift])
    plt.xlabel('Frequency [MHz]')

    frame.axes.get_yaxis().set_visible(False)
    #plt.title(caption.replace('_', ' '))

    plt.savefig(pltdir + 'overlay_' + caption + '.png',
                bbox_inches='tight')
    plt.close()

def plot_scan_with_hist(freq, residuals, prediction=None, label=None, caption=''):
    """
    Show the scatter plot of fluctations vs frequency as well as the distribution
    of the fluctuations. Working from this example:
    http://matplotlib.1069221.n5.nabble.com/how-to-make-scatter-plot-and-
    bar-graphs-in-same-figure-td20364.html
    """
    fig = plt.Figure( (8,8) )
    left, width = 0.1, 0.65 
    bottom, height = 0.1, 0.65 
    bottom_h = left_h = left+width+0.02 

    rect1 = [left, bottom, width, height] 
    rect2 = [left, bottom_h, width, 0.2] 
    rect3 = [left_h, bottom, 0.2, height]
    axScatter = plt.axes(rect1)
    axScatter.scatter(freq, residuals* 1.e24, label=label)

    if prediction is not None:
        axScatter.scatter(freq, prediction * 1.e24,
                          label='predicted signal', color='k')

    axScatter.set_xlabel('Frequency [MHz]')
    axScatter.set_ylabel('Power [yoctooWatts]')
    axScatter.legend(scatterpoints=1)

    nullfmt   = NullFormatter()
    axHisty = plt.axes(rect3)#, sharey=axScatter)

    #formatting
    axHisty.yaxis.set_major_formatter(nullfmt)
    plt.setp(axScatter.get_xticklabels()[-2:], visible=False)
    plt.setp(axScatter.get_xticklines()[-2:], visible=False)
    plt.setp(axHisty.get_yticklabels(), visible=False)
    axHisty.yaxis.set_tick_params(size=0)

    axHisty.hist(residuals * 1.e24, bins=30, orientation='horizontal')

    axHisty.set_ylim( axScatter.get_ylim() )
    for tl in axHisty.get_yticklabels():
        tl.set_visible(False)

    plt.savefig(pltdir + 'scanhist_' + caption + '.png', bbox_inches='tight')
    plt.close()
