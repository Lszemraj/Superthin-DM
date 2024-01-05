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

# Making plots of all five galaxy outputs on one graph
# 3643 was having problems so I removed it from the plot
galaxy_name_array = np.array(["191", "719", #"3643",
                              "4106"])
m_disk_array = np.array([10.580071915555799, 10.205200258802087, #10.58449320302808,
                         10.376780790286702])
m_200_array = np.array([12.248140226879482, 11.842319228037873, #13.78418269561286,
                        12.533385416745602])
log_c_array = np.array([1.006745356947675, 1.1576042848817, #0.6446590844109489,
                        0.9969441832715394])
m_disk_errs_array = np.array([[0.03632496691918696, 0.029021450696498974], [0.26199034515095754, 0.285343323713521],
                             # [ 0.05297980905469579, 0.056903158216202954],
                              [0.17116904040404712, 0.17247754497703305]])
m_200_errs_array = np.array([[0.1741510719374464,0.1754137058092926], [0.3533802066097156, 0.2016428617482653],
                            # [0.11451869275313697, 0.10146411349804829],
                             [0.5441915768588892, 0.8182405503289107]])
log_c_errs_array = np.array([[0.08927741054146487, 0.09719910632486606], [0.09830302518808764, 0.13497510966668402],
                            # [0.03360581549402342, 0.022066310082573426],
                             [0.3368815494850995, 0.17663216767069712]])
colors = ['red', 'orange', #'forestgreen',
          'teal'] #, 'navy']
# plot of m200 vs mstar relation
fig, ax1 = plt.subplots()
m200r = np.arange(9,15,0.1)
ax1.plot(m200r, mass_mhalo(m200r), label='Yang Parametrization')

plt.errorbar(m_200_array, m_disk_array, yerr=m_disk_errs_array.T, xerr=m_200_errs_array.T, fmt = 'o', capsize=2,
             capthick=2, color = "black", ecolor="orange")
ax1.scatter(m_200_array, m_disk_array,marker='*', color = colors)

plt.xlabel("$\log (m200)$")
plt.ylabel("Mstar")
for i, txt in enumerate(galaxy_name_array):
    plt.annotate(txt, (m_200_array[i], m_disk_array[i]))
plt.legend()
plt.title("Plot of m200 vs mstar relation")
plt.show()


# plot of m200-concentration relation
fig2, ax2 = plt.subplots()
ax2.plot(m200r, mass_conc(m200r), label='Yang Parametrization')
#ax2.plot(m_200, log_c,marker='*', ls='none', ms=10, label = f'{galaxy_name} Fit')
plt.errorbar(m_200_array, log_c_array, yerr=log_c_errs_array.T, xerr=m_200_errs_array.T, fmt = 'o', capsize=2,
             capthick=2, color = "black", ecolor="orange")

for i, txt in enumerate(galaxy_name_array):
    plt.annotate(txt, (m_200_array[i], log_c_array[i]))

plt.xlabel("$\log (m200)$")
plt.ylabel("$\log (mass.concentration(m200)$")
plt.legend()
plt.title("Plot of m200-concentration relation")

plt.show()