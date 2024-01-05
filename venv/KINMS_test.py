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

# Working through Kinms tutorial

fit = kinms_fitter('/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/Output_simcube.fits',spatial_trim=[16,48,16,48],spectral_trim=[43,98],linefree_chans=[30,45])

fit.pa_guess = 90
fit.nrings = 5
fit.rms=1e-3 # Jy/beam

# simple fit example
# bestvals, besterrs, acceptedsteps, acceptedsteps_ll, parsfixed = fit.run(method='simple')

# truths=[90,12,10,0,60,25,10,5,  100., 210., 210., 210., 210.] # the true values we would hope to retrieve
# print(Table([fit.labels,bestvals,truths],names=('Quantity', 'Bestfit', 'True value')))


# better fit example using markov chains

fit.niters = 30000
#fit.initial_guesses = bestvals # start from our previous best fit
#fit.show_corner = True # show the corner plots
#bestvals, besterrs, acceptedsteps, acceptedsteps_ll, parsfixed = fit.run(method='mcmc')


# surface brightness profiles

#guesses=[10]
#mymodel=[sb_profs.expdisk(guesses=guesses,minimums=[0],maximums=[20],fixed=[False])]

r=np.arange(0,30,0.1)
#plt.plot(r,sb_profs.eval(mymodel,r,guesses))
#plt.xlabel("Radius")
#plt.ylabel("Brightness")
#plt.show()

#mymodel=[sb_profs.expdisk(guesses=[10],minimums=[0],maximums=[20]),\
        # sb_profs.expdisk(guesses=[1,3],minimums=[0,0],maximums=[100,20],fixed=[False,False])]



# now adding a gaussian
#mymodel=[sb_profs.expdisk(guesses=[10],minimums=[0],maximums=[20]),\
         #sb_profs.gaussian(guesses=[2,10,1],minimums=[0,0,0],maximums=[10,100,20])]

# now adding a cutoff
#mymodel=[sb_profs.expdisk(guesses=[10],minimums=[0],maximums=[20]),\
        # sb_profs.gaussian(guesses=[2,10,1],minimums=[0,0,0],maximums=[10,100,20]),\
         #sb_profs.cutoff(guesses=[0,5],minimums=[0,1],maximums=[10,20])]


#guesses=np.concatenate([i.guess for i in mymodel]) # this line just sets the parameters to the guesses inputed above
#plt.plot(r,sb_profs.eval(mymodel,r,guesses),label='Total')
#plt.plot(r,sb_profs.eval([mymodel[0]],r,guesses[0:1]),label='Comp 1')
#plt.plot(r,sb_profs.eval([mymodel[1]],r,guesses[1:]),label='Comp 2')
#plt.legend(frameon=False)
#plt.xlabel("Radius")
#plt.ylabel("Brightness")
#plt.show()


# velocity profiles

bincentroids=np.arange(0,30,5) # radii of the bin centroids

# initial
#mymodel = [velocity_profs.tilted_rings(bincentroids,guesses=[0,250,260,260,250,230],minimums=np.zeros(bincentroids.size)
                                       #,maximums=np.ones(bincentroids.size)*500.)]

# adding SMBH
#mymodel = [velocity_profs.tilted_rings(bincentroids,guesses=[0,250,260,260,250,230],minimums=np.zeros(bincentroids.size),maximums=np.zeros(bincentroids.size)),\
       #  velocity_profs.keplarian(distance=16.5,guesses=[8],minimums=[5],maximums=[10])]

# replacing tilted rings with arctan
#mymodel = [velocity_profs.arctan(guesses=[240,2],minimums=[0,0],maximums=[500,20]),\
         #velocity_profs.keplarian(distance=16.5,guesses=[8],minimums=[5],maximums=[10])]

# MGE model for star potential

surf=np.array([25551.5, 21118.8, 7436.97, 12016.7, 5862.67, 741.344, 807.669, 212.118])
sigma=np.array([0.226508, 0.661430, 1.30613, 2.17346, 4.76300, 11.3177, 19.2433, 48.5786])
qobs=np.array([0.514866, 0.607566, 0.887495, 0.576108, 0.837162, 0.440516, 0.779643, 0.821153])

mymodel=[velocity_profs.mge_vcirc(surf,sigma,qobs,16.5,guesses=[2,8],minimums=[0.1,7],maximums=[10,10])]


params=np.concatenate([i.guess for i in mymodel])
plt.plot(r,velocity_profs.eval(mymodel,r,params))
#plt.plot(r,velocity_profs.eval([mymodel[0]],r,params[0:-1]),label='Comp 1')
#plt.plot(r,velocity_profs.eval([mymodel[1]],r,params[-1:]),label='Comp 2')
plt.xlabel("Radius")
plt.ylabel("Velocity (km/s)")
plt.show()