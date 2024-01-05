from kinms_fitter import kinms_fitter
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
from astropy.table import Table
from kinms_fitter.sb_profs import sb_profs
from kinms_fitter.velocity_profs import velocity_profs
from kinms_fitter import prior_funcs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString

#plt.step
#find peak flux, then find channels where greater than 50 peak flux
#w50 = 2vmax*sini
#plot error band of W50 -+ deltaV
#decide w20&w50

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]

vel_radio = np.loadtxt(fname = "/Users/Lillie/Desktop/Research Summer 2023/HI_spectra/RFGC4106_HIspec_reduced.txt", skiprows=4, dtype = float, usecols=0)
flux_density = np.loadtxt(fname = "/Users/Lillie/Desktop/Research Summer 2023/HI_spectra/RFGC4106_HIspec_reduced.txt", skiprows=4, dtype = float, usecols=1)

np.set_printoptions(precision = 15)

fig, ax1 = plt.subplots()
max_val = np.max(flux_density)
_thres = 0.2*max_val # width at percentage of peak height
plt.step(vel_radio, flux_density)
plt.axhline(_thres, color = 'red')
plt.xlabel("Velocity (radio)")
plt.ylabel("Flux Density")
plt.legend()
plt.title("H1 Line Test")
plt.show()

thres_array =  np.full((len(vel_radio), 1), _thres)

'''
closest = vel_radio[np.where(flux_density>0.2*np.nanmax(flux_density))]
print(closest)
width = abs(closest[0] - closest[-1]) # width of band (found from first and last intersection)
print("width = ", width) # 458.68111029956617
'''


#finding intersection
first_line = LineString(np.column_stack((vel_radio, flux_density)))
second_line = LineString(np.column_stack((vel_radio, thres_array)))
intersection = first_line.intersection(second_line)
points = [p for p in intersection.geoms]
listarray = []
for pp in points:
    listarray.append([pp.x, pp.y])
nparray = np.array(listarray) # array of point vals, with first element of each point being the vals we need
print(nparray)
width = nparray[0][0] - nparray[-1][0] # width of band
print("width = ", width) #  463.6716592662915
