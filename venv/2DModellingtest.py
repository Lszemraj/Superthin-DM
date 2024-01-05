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


fits_file_path = "/Users/Lillie/PycharmProjects/pythonProject/FitsTester/venv/NGC1387_mom1.fits"


hdul = fits.open(fits_file_path)
hdr = hdul[0].header
data_image = hdul[0].data.T
data_image -= np.nanmedian(data_image)
naxis = hdr['NAXIS'] # should be 3 (number of axes)
naxis1 = hdr['NAXIS1'] # Number of positions along axis 1
naxis2 = hdr['NAXIS2'] # Number of positions along axis 2
crval1 = hdr['CRVAL1'] # value of the physical coordinate identified by CTYPEn at the reference point on axis 1
crval2 = hdr['CRVAL2'] # value of the physical coordinate identified by CTYPEn at the reference point on axis 2
 # value of the physical coordinate identified by CTYPEn at the reference point on axis 3
cdelt1 = hdr['CDELT1']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting
# index,evaluated at the reference point.
cdelt2 = hdr['CDELT2']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting
# index, evaluated at the reference point.
crpix1 = hdr['CRPIX1'] # location along axis n called the reference pixel, or reference point, used in defining the range
# of values for the physical coordinate of axis n.
crpix2 = hdr['CRPIX2']




# printing all file names
# with fits.open(fits_file_path) as hdulist:  # this is like the "hdulist = fits.open('test1.fits')"
#      hdulist.info()
#      for hdu in hdulist:
#           print(repr(hdu.header))


ex_x, ex_y = np.meshgrid(np.arange(- naxis2 /2 , naxis2 / 2,1), np.arange(-naxis1 / 2, naxis1 / 2,1))


#truth = np.array([crval3, 0.5, 5, 270, crpix1,crpix2])
#ex_mom1 = moment1(truth, ex_x, ex_y)


mcmc = gastimator(moment1,ex_x*cdelt1, ex_y*cdelt2)
mcmc.labels = np.array(['Vmax','sin_I', 'rturn', 'PA', 'xc', 'yc'])
mcmc.guesses = np.array([160, 0.5, 1, 45, 4, 4])
mcmc.min = np.array([50.,0., 0., 0., -5., -5.])
mcmc.max = np.array([250.,1., 2., 359., 5., 5.])
mcmc.fixed = np.array([False, True, False, False, False, False])

mcmc.precision = np.array([10.,5., 1.,5.,2.,2.])
nsamples = 100000
error = 10.
outputvalue, outputll= mcmc.run(data_image,error,nsamples,nchains=7,plot=False)
bestfit_model= moment1(mcmc.guesses, ex_x, ex_y)
# print("model values",np.median(bestfit_model))
figure = corner_plot.corner_plot(outputvalue[~mcmc.fixed,:].T,like=outputll,labels=mcmc.labels[~mcmc.fixed],quantiles=[0.16, 0.5, 0.84])
# fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
# a = ax1.contourf(ex_x, ex_y, data_image,levels= np.arange( -150, 150,10))
# b = ax2.contourf(ex_x, ex_y, bestfit_model,levels= np.arange( -150, 150,10))

plt.show()