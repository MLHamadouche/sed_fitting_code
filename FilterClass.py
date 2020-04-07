import numpy as np
import matplotlib.pyplot as plt
import pickle

#give array/list of filter paths

class Filters:

    def __init__(self, filter_list, wavelengths):
        self.filter_list = filter_list
        self.wavelengths = wavelengths
        #self.fluxes  = fluxes

        self.load_filter_files()
        self.calc_eff_wavs()
        #self.filter_interp()


    def load_filter_files(self):
        filter_curves = []

        for filter in filter_list:
            filter_curves.append(np.loadtxt("/Users/PhDStuff/sed_fitting_code/filters/"+str(filter)))

        self.filter_curves = filter_curves
            #print(self.filter_curves)
    def calc_eff_wavs(self):
        eff_wavs = []

        for f in self.filter_curves:
            flux_filter = f[:,1]/np.max(f[:,1])
            wav_filter = f[:,0]
            eff_wavs.append(np.sum(wav_filter*flux_filter)/np.sum(flux_filter))

        self.eff_wavs = eff_wavs
        #print(self.eff_wavs)

    #def filter_interp(self, redshift, wavelengths):
        #new_filter_curves = np.zeros(len(self.filter_list))

        #for f in self.filter_curves:
        #    flux_filter = f[:,1]/np.max(f[:,1])
            #wav_filter = f[:,0]
            #filter_interpol = np.interp(wavelengths*(1 + redshift), wav_filter, flux_filter, left=0, right=0) # shift Wavelengths by factor of redshift to the left to get lambda_obs
        #    new_filter_curves=filter_interpol

        #self.new_filter_curves = new_filter_curves


file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
model = flux_grid[4][150, :]

filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]
#print(filter_list)

waves = np.array(waves)
redshift = 0.1
fluxes = model
filters = Filters(filter_list, waves)
