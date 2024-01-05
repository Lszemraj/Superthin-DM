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
cdelt1 = hdr['CDELT1']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting index,
# evaluated at the reference point.
cdelt2 = hdr['CDELT2']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting index,
# evaluated at the reference point.
 # rate of change of the physical coordinate along axis n per unit change in the counting index,
# evaluated at the reference point.
crpix1 = hdr['CRPIX1'] # location along axis n called the reference pixel, or reference point, used in defining the range
# of values for the physical coordinate of axis n.
crpix2 = hdr['CRPIX2']


# approx vals from fitting
Vmax = 133.4756114414613
sin_I = 0.5
rturn = 0.12993552360122673
PA = 334.74876573780614
xc = -0.3031411141469313
yc = 0.23015565447990669
dist=19.9#Mpc
#array_axis1 = ((np.arange(0, naxis1) - (crpix1 - 1)) * cdelt1) #+ crval1
#array_axis2 = ((np.arange(0, naxis2) - (crpix2 - 1)) * cdelt2) #+ crval2
#array_axis1,array_axis2 = np.arange(- naxis2 /2 , naxis2 / 2,1), np.arange(-naxis1 / 2, naxis1 / 2,1)

r = np.arange(0,110) #np.sqrt((array_axis1-xc)**2 + (array_axis2-yc)**2)
print(r.min(),r.max())
Vr = 2 * (Vmax / math.pi) * np.arctan(r / rturn)


# plot velocity curve graph
fig, ax = plt.subplots()
ax.plot(r*np.abs(cdelt1)*4.84*dist, Vr)
ax.set_title("Velocity Curve Graph")
ax.set_xlabel("Radius in pc")
ax.set_ylabel("Velocity in km/s")
plt.show()


# mass enclosed equation
def enclosed_mass_function(r):
    Vr = 2 * (Vmax / math.pi) * np.arctan(r / rturn)
    G = 4.4301e-3 # in units of Msolar pc^-1 km/s^-2
    return (r*(Vr**2))/G

examined_r = 60
print(f"Enclosed mass at radius {examined_r} pc in M_solar:", enclosed_mass_function(examined_r))


fig, ax = plt.subplots()
ax.semilogy(r*np.abs(cdelt1)*4.84*dist, enclosed_mass_function(r*np.abs(cdelt1)*4.84*dist))

plt.show()


