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

# 191
distance = 82.03 #in Mpc
r = np.arange(0,160,2)
vel = 437.51928131017485 / 2
galaxy_name = "RFGC191"

mymodel = [velocity_profs.exponential_disc(distance,guesses=[10.580071915555799,3.592533293026054],minimums=[6,1],maximums=[12.5,20])
    , velocity_profs.nfw(distance,guesses=[12.248140226879482,1.006745356947675],minimums=[7,0],maximums=[16,2])]



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



