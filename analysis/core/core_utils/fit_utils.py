import numpy as np

import scipy.optimize as optimize


def chi2(p, x, y, fp, err):
    """Compute chi squared."""
    chi2 = 0.
    
    for i in range(len(x)):
        model = fp(x[i])
        chi2 += (model - y[i])**2 / err[i]**2
    
    return chi2
    
def simplex_fit(x, y, err=None, deg=5):
    """
    Use the simplex (Nealder-Mead) algorithm to fit the data.
    Simplex is more robust than Levenberg-Marquardt, but
    does not return uncertainties for the estimate.
    Assume the parametric function is a deg order polynomial.
    https://glowingpython.blogspot.com/2011/05/curve-fitting
    -using-fmin.html
    """
    #parametric function of degree deg with c a list of the
    #parameters
    fp = lambda c, x: sum((x**i) * coeff for i, coeff
                                  in enumerate(c))

    #error function to minimize. If uncertainties are
    #passed as input, use chi squared; otherwise, use
    #a simple sum of errors.
#    if err is not None:
#        e = lambda p, x, y, err: chi2(p, x, y, fp, err)
    if True:
        e = lambda p, x, y: abs(fp(p, x) - y).sum()
    
    #initial parameters
    guess = [0] * (deg + 1)
    guess[0] = y.mean() / 2.
    guess[1] = y.mean() / 2. * 1. / x.mean()

    #fit the data
    p = optimize.fmin(e, guess, args=(x,y))

    return fp(p, x), p
    
def  poly_fit(x, y, deg=5):
    """
    Use numpy's polyfit. Does not return uncertainty
    for the estimate.
    """
    fit_params = np.polyfit(x, y, deg)
    y_new = np.polyval(fit_params, x)

    return y_new, fit_params

def leastsq_fit(x, y, deg=5):
    """
    Use scipy's curve_fit, a wrapper for it's leastsquares
    algorithm, which is a modified Levenberg-Marquardt alg.
    """
    #parametric function of degree deg with c a list of the
    #parameters. Note the order is x then c!
    def fp(c, x):
        output = 0
        for i, coeff in enumerate(c):
            output += coeff * x**i
        return output

    #error function whose square will be minimized by leastsq
    e = lambda p, x, y: y - fp(p, x)

    #initial parameters
    guess = [0] * (deg + 1)
    guess[0] = y.mean() / 2.
    guess[1] = y.mean() / 2. * 1. / x.mean()

    #(TODO) should I scale the data to 1? I've heard curve_fit
    #has an easier time if data is not extremely small or large.

    fit_params, cov = optimize.leastsq(e, guess, args=(x,y))
    
    return fp(fit_params, x), fit_params

if __name__=="__main__":
    import matplotlib.pyplot as plt

    order = 2
    #test out a known polynomial of degree = order
    x = np.linspace(6.e8, 6.001e8, 30)
    fp = lambda c, x: sum((x**i) * coeff for i, coeff
                          in enumerate(c))
    scale = 1.e-5
    real_p = np.random.rand(3)
    y = fp(real_p, x) + np.random.normal(0, 0.05, 30)

    plt.plot(x, y*scale, 'bo')
    for func in [simplex_fit, poly_fit]:
        result, params = func(x, y*scale, deg=order)
        plt.plot(x, result, label=str(func).split('_')[0])
        print params

    plt.legend()
    plt.show()
