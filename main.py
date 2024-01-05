import astropy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

fits_file_path = (
    "/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/NGC6753_mom0.fits"
)
hdul = fits.open(fits_file_path)
hdr = hdul[0].header
data_image = hdul[0].data

# printing all file names
with fits.open(
    fits_file_path
) as hdulist:  # this is like the "hdulist = fits.open('test1.fits')"
    hdulist.info()
    for hdu in hdulist:
        print(repr(hdu.header))


# accessing elements of the header
naxis = hdr["NAXIS"]  # should be 2 (number of axes)
naxis1 = hdr["NAXIS1"]  # Number of positions along axis 1
naxis2 = hdr["NAXIS2"]  # Number of positions along axis 2
crval1 = hdr["CRVAL1"]
crval2 = hdr["CRVAL2"]
cdelt1 = hdr["CDELT1"]
cdelt2 = hdr["CDELT2"]
crpix1 = hdr["CRPIX1"]
crpix2 = hdr["CRPIX2"]

# coordinate thingy
array_axis1 = ((np.arange(0, naxis1) - (crpix1 - 1)) * cdelt1) + crval1
array_axis2 = ((np.arange(0, naxis2) - (crpix2 - 1)) * cdelt2) + crval2


# plot fits file
fig, ax = plt.subplots(1, 1, figsize=(7, 7))
ax.contourf(
    array_axis1,
    array_axis2,
    data_image,
    levels=np.linspace(0, 1, 10) * np.nanmax(data_image),
)
plt.show()
