import numpy as np
from astropy.io import fits
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


file = "/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/hst_10217_09_acs_wfc_f850lp_drz.fits"
f = fits.open(file)
data = f[1].data
hdr = f[1].header


# puts image in AB magnitudes
instr_ABMAG_ZPT =-2.5*np.log10(hdr['PHOTFLAM'])-21.10-(5*np.log10(hdr['PHOTPLAM']))+18.6921
data_ABmag =- 2.5*np.log10(data)+instr_ABMAG_ZPT


naxis = hdr['NAXIS'] # should be 3 (number of axes)
naxis1 = hdr['NAXIS1'] # Number of positions along axis 1
naxis2 = hdr['NAXIS2'] # Number of positions along axis 2
crval1 = hdr['CRVAL1'] # value of the physical coordinate identified by CTYPEn at the reference point on axis 1
crval2 = hdr['CRVAL2'] # value of the physical coordinate identified by CTYPEn at the reference point on axis 2
cdelt1 = hdr['CDELT1']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting
# index,evaluated at the reference point.
cdelt2 = hdr['CDELT2']*3600 # rate of change of the physical coordinate along axis n per unit change in the counting
# index,evaluated at the reference point.
crpix1 = hdr['CRPIX1'] # location along axis n called the reference pixel, or reference point, used in defining the
# range of values for the physical coordinate of axis n.
crpix2 = hdr['CRPIX2']


dist = 19.9 # in pc
xc = -0.3031411141469313
yc = 0.23015565447990669

ex_x, ex_y = np.arange(- naxis2 /2 , naxis2 / 2,1), np.arange(-naxis1 / 2, naxis1 / 2,1)
r = np.arange(0,4969) #r = np.sqrt((ex_x-xc)**2 + (ex_y-yc)**2)

Vmax = 133.4756114414613
sin_I = 0.5
rturn = 0.12993552360122673
PA = 334.74876573780614
xc = -0.3031411141469313
yc = 0.23015565447990669
dist=19.9#Mpc

ab_mag_at_r = np.sqrt((data_ABmag[0])**2 + (data_ABmag[1])**2)

def calc_light(r, data, center_x, center_y):
    sizes = data.shape
    r_squared = r**2
    sum = 0
    for x in range(sizes[0]):
        for y in range(sizes[1]):
            dist_x = x - center_x
            dist_y = y - center_y
            new_r = dist_x**2 + dist_y**2
            if new_r <= r_squared:
                sum += data[x][y]
    return sum



# mags_summed = []
# for i in range(0, 400):
 # mags_summed.append(calc_light(i, data_ABmag, 2813, 2540))

def enclosed_mass_function(r):
    Vr = 2 * (Vmax / math.pi) * np.arctan(r / rturn)
    G = 4.4301e-3 # in units of Msolar pc^-1 km/s^-2
    return (r*(Vr**2))/G



#data_ABmag_copy = data_ABmag[2300:2700,2600:3000]

#fig, ax = plt.subplots()
#ax.plot(r*np.abs(cdelt1)*4.84*dist, ab_mag_at_r)
#plt.pcolormesh(data_ABmag)
#plt.show()

f=fits.open(file)
data=f[1].data[2813-200:2813+200,2540-200:2540+200]
hdr=f[1].header
cellsize=hdr['CD2_2']*3600.
s=data.shape
x=(np.arange(0,s[0])-200)*cellsize
y=(np.arange(0,s[1])-200)*cellsize



ygrid,xgrid=np.meshgrid(x,y)

rsqr=(ygrid**2 + xgrid**2)

rflat=np.arange(0,10.35,0.1)
massr= np.zeros_like(rflat)

for i,rval in enumerate(rflat):
    massr[i]=np.nansum(data[rsqr<rval**2])


plt.semilogy(rflat*4.84*dist, (massr/massr.max()) * 1.52122771e+10 * 0.25) # in units of M sun
plt.semilogy(rflat*4.84*dist, enclosed_mass_function(rflat*4.84*dist))
plt.xlabel("R [pc]")
plt.ylabel("Mass Enclosed [M$_{\odot}$]")
plt.show()
