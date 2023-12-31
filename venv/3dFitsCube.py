import astropy
from astropy.io import fits
import numpy as np
from numpy import zeros, newaxis
import matplotlib.pyplot as plt



fits_file_path = "/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/NGC1387_cube.fits"
hdul = fits.open(fits_file_path)
hdr = hdul[0].header
data_image = hdul[0].data.T

# printing all file names
with fits.open(fits_file_path) as hdulist:  # this is like the "hdulist = fits.open('test1.fits')"
     hdulist.info()
     for hdu in hdulist:
          print(repr(hdu.header))


# accessing elements of the header
naxis = hdr['NAXIS'] # should be 3 (number of axes)
naxis1 = hdr['NAXIS1'] # Number of positions along axis 1
naxis2 = hdr['NAXIS2'] # Number of positions along axis 2
naxis3 = hdr['NAXIS3']
crval1 = hdr['CRVAL1']
crval2 = hdr['CRVAL2']
crval3 = (hdr['CRVAL3'] / 1000)
cdelt1 = hdr['CDELT1']
cdelt2 = hdr['CDELT2']
cdelt3 = (hdr['CDELT3'] / 1000)
crpix1 = hdr['CRPIX1']
crpix2 = hdr['CRPIX2']
crpix3 = hdr['CRPIX3']

# correction
array_axis1 = ((np.arange(0, naxis1) - (crpix1 - 1)) * cdelt1) + crval1
array_axis2 = ((np.arange(0, naxis2) - (crpix2 - 1)) * cdelt2) + crval2
array_axisV = ((np.arange(0, naxis3) - (crpix3 - 1)) * cdelt3) + crval3

data_image[np.isfinite(data_image) == False] = 0
data_image[data_image < 5*np.nanstd(data_image[:,:,0:5])] = 0
momentZero = np.nansum(data_image, axis = 2)
momentZeroNew = momentZero[:, :, newaxis]

momentOne = (data_image * array_axisV) / momentZeroNew
momentOnesummed = np.nansum(momentOne, axis = 2)
#print(momentOnesummed)
#print(momentOnesummed.shape)

print(np.arange(-200,200,cdelt3))
# momentZero plot
fig,ax= plt.subplots(1,1,figsize=(7,7))
v_min = np.nanmin(momentOnesummed)
v_max = np.nanmax(momentOnesummed)
v_med = np.nanmedian(momentOnesummed[momentOnesummed != 0])
print("v_med: ", v_med)
print("array_v :", array_axisV )
a = ax.contourf(array_axis1, array_axis2, momentOnesummed,levels= np.arange(v_med - 200,v_med + 200,cdelt3))
plt.colorbar(a)
plt.show()