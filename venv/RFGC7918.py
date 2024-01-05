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

class Yang_Model():
    """
    Parameterized SMHM relation from Yang et al. 2003.

    Default parameters are from Moster et al. 2014 table 1.  They have been altered to be
    in h=1 units.

    Parameters
    ==========
    norm: float, optional
        characteristic stellar-to-halo mass ratio

    M1: float, optional
        log characteristic halo mass

    beta: float, optional
        low mass slope

    gamma: float, optional
        high mass slope

    """

    def __init__(self, norm=0.0203, M1=11.741, beta=1.057, gamma=0.556):
        self.norm = norm
        self.M1 = M1
        self.beta = beta
        self.gamma = gamma

    def __call__(self, Mh):
        """
        Parameters
        ==========
        Mh: array_like
            log halo mass

        Returns
        =======
        Mstar: array_like
            log stellar mass

        """
        Mh = 10.0 ** Mh
        M1 = self.M1
        norm = self.norm
        beta = self.beta
        gamma = self.gamma
        M1 = 10.0 ** M1
        mstar = 2.0 * norm * ((Mh / M1) ** (-1.0 * beta) + (Mh / M1) ** gamma) ** (-1.0) * Mh
        return np.log10(mstar)


mass_mhalo = Yang_Model()


def mass_conc(m200):
    return 1.025 - 0.097 * (m200 - 12)



# load in galaxy
cube  = '/Users/Lillie/Desktop/Research Summer 2023/UltraThin.Galaxies/RFGC719.CO.cube.fits'
distance = 126.20142857142858


#plotter = pymakeplots(cube= cube)
#plotter.make_all(fits = True, pdf = True)


# do fitting process
fit= kinms_fitter(cube, spatial_trim = [131-60,131+60,131-60,131+60],spectral_trim = [210,270])
fit.pa_guess= 530.
fit.inc_guess= 90.
fit.rms=1e-3 # Jy/beam
fit.vsys_guess = 8500+85
fit.totflux_guess = 4.2
massprior = priors.gaussian(10.19,0.3)

fit.sb_profile = [sb_profs.expdisk(guesses=[4],minimums=[0.1],maximums=[10])]

fit.vel_profile=[velocity_profs.exponential_disc(distance,guesses=[10.19,4.3],minimums=[4,1],maximums=[14.5,20],priors=(massprior.eval,None))
    , velocity_profs.nfw(distance,guesses=[11.76,.9],minimums=[5,0],maximums=[14,2])]
fit.niters = 10000
fit.pdf = True
#bestvals, besterrs, outputvalue, outputll,_ = fit.run(method='both' )#, justplot = 'True')


# properties from best fit
m_disk =  10.205200258802087
m_200 = 11.571104316881046
m_disk_error_up = 0.26199034515095754
m_disk_error_down = 0.285343323713521
m_200_error_up = 0.3533802066097156
m_200_error_down = 0.2016428617482653
log_c = 1.1576042848817
log_c_error_up = 0.09830302518808764
log_c_error_down = 0.13497510966668402

galaxy_name = 'RFGC719'

# plot of m200 vs mstar relation
fig, ax1 = plt.subplots()
m200r = np.arange(9,15,0.1)
ax1.plot(m200r, mass_mhalo(m200r), label='Yang Parametrization')
#ax1.plot(m_200, m_disk,marker='*', ls='none', ms=10, label = f'{galaxy_name} Fit')
plt.errorbar(m_200, m_disk, yerr=np.array([[m_disk_error_up, m_disk_error_down]]).T, xerr=np.array([[m_200_error_up, m_200_error_down]]).T)
plt.xlabel("$\log (m200)$")
plt.ylabel("Mstar")
plt.legend()
plt.title("Plot of m200 vs mstar relation")
plt.show()

# plot of m200-concentration relation
fig2, ax2 = plt.subplots()
ax2.plot(m200r, mass_conc(m200r), label='Yang Parametrization')
#ax2.plot(m_200, log_c,marker='*', ls='none', ms=10, label = f'{galaxy_name} Fit')
plt.errorbar(m_200, log_c, yerr=np.array([[log_c_error_up, log_c_error_down]]).T, xerr=np.array([[m_200_error_up, m_200_error_down]]).T)
plt.xlabel("$\log (m200)$")
plt.ylabel("$\log (mass.concentration(m200)$")
plt.legend()
plt.title("Plot of m200-concentration relation")

plt.show()