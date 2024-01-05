import math
from pymakeplots import pymakeplots
import astropy
from astropy.io import fits
import numpy as np
from numpy import zeros, newaxis
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from gastimator import gastimator
from gastimator import corner_plot

# write function to spit out moment one (
def moment1(vals, x, y):
    Vmax = vals[0]
    sin_I = vals[1]
    rturn = vals[2]
    PA = vals[3]
    xc = vals[4]
    yc = vals[5]
    #sin_I=0.5
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    Vr = 2 * (Vmax / math.pi) * np.arctan(r / rturn)
    theta =  np.arctan2(x-xc, y-yc) + np.deg2rad(PA)
    mom1 =  Vr * sin_I * np.cos(theta)
    return mom1

#function for chi squared
def chi_squared(observed, expected):
    x = (observed - expected)**2 / expected
    x = np.nansum(x)
    return x

#minimize with curve fit
#scipy.curve_fit()


ex_x, ex_y = np.meshgrid(np.arange(-25, 25,1), np.arange(-25, 25,1))
truth = np.array([200, 0.5, 5, 270,0,0])
ex_mom1 = moment1(truth, ex_x, ex_y)
#plt.imshow(ex_mom1, origin="lower")
#plt.show()

#using scipy optimize curve fit
#optimize.curve_fit(moment1,(ex_x, ex_y), ex_mom1)


# using GAStimator

mcmc = gastimator(moment1,ex_x, ex_y)
mcmc.labels = np.array(['Vmax','sin_I', 'rturn', 'PA', 'xc', 'yc'])
mcmc.guesses = np.array([215, 0.5, 7, 243, 0, 0])
mcmc.min = np.array([190.,0., 4., 0., -5., -5.])
mcmc.max = np.array([220.,1., 8., 360., 5., 5.])
mcmc.fixed = np.array([False, True, False, False, False, False])

mcmc.precision = np.array([10.,5., 1.,5.,2.,2.])
nsamples = 100000
data = ex_mom1
error = 10.
outputvalue, outputll= mcmc.run(data,error,nsamples,nchains=7,plot=False)

figure = corner_plot.corner_plot(outputvalue[~mcmc.fixed,:].T,like=outputll,labels=mcmc.labels[~mcmc.fixed],quantiles=[0.16, 0.5, 0.84], truths=truth[~mcmc.fixed])
plt.show()