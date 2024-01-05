from pymakeplots import pymakeplots
import astropy
from astropy.io import fits
import numpy as np
from numpy import zeros, newaxis
import matplotlib.pyplot as plt

file= '/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/RFGC191.CO.cube.fits'
# file hosted on the web in this case, but pointing the code to a cube on your local system will work.

plotter = pymakeplots(cube= file)
plotter.make_all(fits = True, pdf = True)
