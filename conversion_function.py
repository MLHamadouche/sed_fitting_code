import numpy as np

# convert microJanskys to erg s^-1 cm^-2 A^-1
# bc03 models in erg s^-1 cm^-2 Angstrom^-1

#1 Jy  = 10^-23 erg s^-1 cm^-2 Hz^-1

def conversion_func(fluxes, eff_wavs):
    mu = (10**-6) #cm
    c = 2.99792458*(10**10) #cm/s^2
    fluxes = fluxes*(10**-23)*mu
    fluxes = (fluxes*2.99792458*10**18)/(eff_wavs**2)

    return fluxes
