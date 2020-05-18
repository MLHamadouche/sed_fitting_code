import numpy as np
import astropy
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simps
import cosmology

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
Mpc_m = 3.086*10**22 #Mpc in m
Mpc_cm = Mpc_m*10**2 #in centimetres
#filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]
"""
def load_filter_files(filter_list):
    filter_curves = []

    for filter in filter_list:
        filter_curves.append(np.loadtxt("/Users/PhDStuff/sed_fitting_code/filters/"+str(filter)))

    return filter_curves
"""

def load_filter_files(filter_list):
    filter_curves = []

    for filter in filter_list:
        filter_curves.append(np.loadtxt("/Users/massissiliahamadouche/Downloads/massi_zphot_test/ECDFS_filters/"+str(filter)))
    return filter_curves
#print(self.filter_curves)
def calc_eff_wavs(filter_curves):
    eff_wavs = []

    for f in filter_curves:
        flux_filter = f[:,1]/np.max(f[:,1])
        wav_filter = f[:,0]
        eff_wavs.append(np.sum(wav_filter*flux_filter)/np.sum(flux_filter))

    return eff_wavs

class Photometry:

    def __init__(self, waves, filter_curves, eff_wavs, redshift, fluxes):
        self.waves = waves
        self.filter_curves= filter_curves
        self.redshift = redshift
        self.fluxes = fluxes
        self.eff_wavs = eff_wavs
        #self.load_filter_files()
        #self.calc_eff_wavs()
        self.cosmos()
        self.filter_interp()
        self.photometry()


    def filter_interp(self):
        new_filter_curves = []
        for f in self.filter_curves:
            flux_filter = f[:,1]/np.max(f[:,1])
            wav_filter = f[:,0]
            filter_interpolation = np.interp(self.waves*(1 + self.redshift), wav_filter, flux_filter, left=0, right=0) # shift Wavelengths by factor of redshift to the left to get lambda_obs
            new_filter_curves.append(filter_interpolation)

        self.new_filter_curves = new_filter_curves


    def cosmos(self):

        self.fluxes*=(3.826*10**33) #for conversion to ergs/s/cm^2/Angstrom
        l_dist=cosmo.luminosity_distance(self.redshift).value*Mpc_cm #converting Mpc to cm
        lum_area = 4 *np.pi*((l_dist)**2)
        f_lambda = self.fluxes/(lum_area)
        f_lambda/=(1+self.redshift)

        self.converted_fluxes = f_lambda

    def photometry(self):


        new_fluxes = np.zeros((len(self.new_filter_curves)))

        for m in range(len(self.new_filter_curves)):
            f_lam_x_filters = self.converted_fluxes*self.new_filter_curves[m]*self.waves
            new_fluxes[m] = (np.trapz(f_lam_x_filters, x=self.waves))/(np.trapz(self.new_filter_curves[m]*self.waves, x=self.waves))
        self.model_photometry = new_fluxes

        return self.model_photometry
