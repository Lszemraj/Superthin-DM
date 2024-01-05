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
from gastimator import priors

cube  = '/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/RFGC3643.CO10.fits'
mom0 = fits.open("/Users/Lillie/PycharmProjects/pythonProject/FitsTester/venv/RFGC3643_mom0.fits")[0].data
mom1 = "/Users/Lillie/PycharmProjects/pythonProject/FitsTester/venv/RFGC3643_mom2.fits"

distance = 118.72285714285715 # in Mpc

fit=kinms_fitter(cube, spatial_trim = [129-100,129+100,129-100,129+100],spectral_trim = [80,160])
fit.pa_guess= 75.
fit.inc_guess= 90.
fit.rms=1e-3 # Jy/beam
fit.vsys_guess=8030+85
fit.totflux_guess=9.5
#fit.m200_guess = 15
#fit.concentration_guess = 2
#fit.mdisk_guess = 10
#fit.rexpscale_guess = 8

massprior = priors.gaussian(10.58,0.3)
fit.sb_profile = [sb_profs.expdisk(guesses=[7],minimums=[0.1],maximums=[15])]
#fit.vel_profile = [velocity_profs.exponential_disc(distance = 19.9, guesses=np.array([15,10]),minimums=np.array([8,0.01]),maximums=[200, 70]),
                  # velocity_profs.nfw(distance = 19.9, guesses=np.array([50,4]),minimums=np.array([4,0.01]),maximums=[200,60])]
#fit.vel_profile=[velocity_profs.exponential_disc(distance,guesses=[10,3.6],minimums=[9,1],maximums=[12.5,30])
   # , velocity_profs.nfw(distance,guesses=[10.5,0.90],minimums=[9,0],maximums=[15,2])]

fit.vel_profile=[velocity_profs.exponential_disc(distance,guesses=[10.58,6.3],minimums=[9,1],maximums=[12.5,30], priors=(massprior.eval,None))
    , velocity_profs.nfw(distance,guesses=[12.66,0.98],minimums=[9,0],maximums=[14,2])]

fit.niters = 10000
fit.pdf = True
bestvals, besterrs, outputvalue, outputll,_ = fit.run(method='both') #, justplot= "True")



#print(Table([fit.labels,bestvals,besterrs],names=('Quantity', 'Bestfit','1-sigma error')))



#print(mom0.shape)
#plt.imshow(mom0)
#plt.show()