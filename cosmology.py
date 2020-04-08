import numpy as np
import astropy
from astropy.cosmology import FlatLambdaCDM


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
Mpc_m = 3.086*10**22 #Mpc in m
Mpc_cm = Mpc_m*10**2 #in centimetres

def cosmos(fluxes, redshift):
        fluxes*=(3.826*10**33) #for conversion to ergs/s/cm^2/Angstrom
        l_dist=cosmo.luminosity_distance(redshift).value*Mpc_cm #converting Mpc to cm
        lum_area = 4 *np.pi*((l_dist)**2)
        f_lambda = fluxes/(lum_area)

        return f_lambda
