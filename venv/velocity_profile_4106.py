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
from gastimator import priors

# 4106
distance = 103.47000000000001 #in Mpc
r = np.arange(0,150,2)
vel =  458.68111029956617 / 2
massprior = priors.gaussian(10.53,0.3)
galaxy_name = "RFGC4106"

mymodel = [velocity_profs.exponential_disc(distance,guesses=[10.376780790286702,3.8197009704115525],minimums=[4,1],maximums=[12.5,20],
                                                 priors=(massprior,None))
    , velocity_profs.nfw(distance,guesses=[12.533385416745602,0.9969441832715394],minimums=[5,0],maximums=[14,2])]



params=np.concatenate([i.guess for i in mymodel])
#plt.plot(r,velocity_profs.eval(mymodel,r,params),label="Total")
plt.plot(r*distance*4.84e-3,velocity_profs.eval([mymodel[0]],r,params[0:2]),label='Exponential Disc')
plt.plot(r*distance*4.84e-3,velocity_profs.eval([mymodel[1]],r,params[2:]),label='NFW')
plt.plot(r*distance*4.84e-3,velocity_profs.eval(mymodel,r,params),label='Total')
plt.axhline(y = vel, color = 'r', label = 'From HI')
plt.legend(frameon=False)
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.title(f"Velocity Profile of {galaxy_name}")
plt.show()
