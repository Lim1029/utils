# general gaussian equation, from https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian1D.html  

from scipy.optimize import curve_fit
import numpy as np

def gaussian(v, amplitude, vlsr, width):
    sigma = width / (2*np.sqrt(2*np.log(2)))
    return amplitude * np.exp(-(v-vlsr)**2/(2*sigma**2))

def single_gaussian_fit(x,y,a0,v0,w0):
    try:
        popt, pcov = curve_fit(gaussian, x, y,
                               p0=[a0,v0,w0], # amplitude, vlsr, width
                               )
        a_fit, v_fit, w_fit = popt 
        area = a_fit*w_fit*1.0645
        return a_fit, v_fit, w_fit, area
    except Exception as e:
        print(f"Fitting Failed: {e}")

def single_gaussian_fit_autoinit(x,y):
    # estimate the initial parameters
    # v0 = moment1
    # w0 = moment2
    # a0 = the flux at moment1
    dx = np.abs(x[1]-x[0])
    mom0 = np.nansum(y*dx)
    mom1 = np.nansum(x*y*dx)/mom0
    # y_masked = y>np.std(y[:20]) #simple 1 sigma masking
    y_positive = np.where(y>0,y,-y)
    mom2 = np.sqrt(np.nansum(y_positive*(x-mom1)**2*dx)/mom0)
    # breakpoint()
    # a0 = y[np.where(x>mom1)][0]
    a0 = np.nanmax(y)
    v0, w0 = mom1, mom2/(2*np.sqrt(2*np.log(2)))
    print(a0,v0,w0)
    try:
        popt, pcov = curve_fit(gaussian, x, y,
                               p0=[a0,v0,w0], # amplitude, vlsr, width
                               bounds=(
                                    [ 0,      -np.inf, 0.0],
                                    [ np.inf,  np.inf, np.inf])
                               )
        a_fit, v_fit, w_fit = popt 
        area = a_fit*w_fit*1.0645
        return a_fit, v_fit, w_fit, area
    except Exception as e:
        print(f"Fitting Failed: {e}")
