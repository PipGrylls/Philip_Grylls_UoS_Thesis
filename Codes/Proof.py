"""This is the original proof of concept for unevolved 
SHMF creating the theoritical upper limit"""

import numpy as np
import matplotlib.pyplot as plt

import SEM # Module containing semi-empirical routines

#Subhalomass function parameters macc/M0
Unevolved = {
'gamma' : 0.22,
'alpha' : -0.91,
'beta' : 6,
'omega' : 3,
'a' : 1}

#Returns the Unevolved SHMF from Jiang, van den Bosch.
#Units are Mvir h-1
def dn_dlnX(Parameters, X):
    """
    Caculates subahlo mass funtions
    Args:
        Parameters: Dictonary containg 'gamma', 'alpha', 'beta', 'omega', 'a'
        X: m/M arrays desired subhalo/parenthalo
    Returns:
        dn_dlogX_arr: Numberdensitys per dex #N dex-1
    """
    Part1 = Parameters['gamma'] * np.power(Parameters['a'] * X, 
                                           Parameters['alpha'])
    Part2 = np.exp(-Parameters['beta'] *
                   np.power(Parameters['a']*X, Parameters['omega']))
    dn_dlnX_arr = Part1*Part2
    dn_dlogX_arr = dn_dlnX_arr*2.30
    return dn_dlogX_arr # N dex-1

Binwidth_SDSS = 0.01
# Range of host halo massed to investigate linear will require weighting later
CentralHaloMass = np.arange(12+ np.log10(h), 
                            15+ np.log10(h), 
                            Binwidth_SDSS) #Mvir h-1

# Range of satilite masses to investigate 
# starting at 10^11 as we are looking for satilites above 10^10
SatBin = 0.1
SatHaloMass = np.arange(11+ np.log10(h), 15+ np.log10(h), SatBin)

# Makes m/M as required by our Jing et al
m_M = np.array([SatHaloMass - i for i in CentralHaloMass])

# Unevolved SHMF

# Runs the model from Jing
Out = dn_dlnX(Unevolved, np.power(10, m_M))
# Weight the output to the HMF(central)
Out_Weighted = np.array([thing*HMF_fun(CentralHaloMass[i]- np.log10(h)) 
                         for i, thing in enumerate(Out)])

# Abundance matching
StellarX = SEM.DarkMatterToStellarMass(SatHaloMass- np.log10(h), 
                                       0, Paramaters, ScatterOn = False)
# Masscuts
StellarX_10 = StellarX[StellarX > SatiliteMassCut]
# Integrates
Integrals = np.array([trapz(thing[StellarX > SatiliteMassCut], StellarX_10) 
                      for thing in Out_Weighted])

# f=sum(Bin)/(sum(population)*binwidth)
AnalyticModel = Integrals/(np.sum(Integrals)*Binwidth_SDSS)

# Plotting
plt.plot(CentralHaloMass, AnalyticModel)
plt.savefig("./Figures/AnalyticMax.png")
plt.clf()

"""Proof of concept over"""