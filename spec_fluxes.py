import numpy as np
import matplotlib.pyplot as plt
import time
import astropy.io.fits as fits
import os.path

#if os.path.exists("/bc03/Miles_Atlas/Kroupa_IMF/"):
#    path = os.path.dirname(os.path.realpath(__file__))
path = "bc03/Miles_Atlas/Kroupa_IMF/"
gen_file = str(path) + "bc2003_hr_xmiless_m62_kroup_ssp.ised_ASCII"


age_bc03 = np.loadtxt(gen_file, max_rows = 1, dtype=float)

ages = np.array(age_bc03[1:222])


lambda_bc03 = np.genfromtxt(gen_file, skip_header=6, skip_footer=233, dtype=float)
lambda_ = np.array(lambda_bc03)

wavelengths = lambda_[1:]

def spectrum():
    models = ["m22","m32", "m42","m52", "m62", "m72", "m82"]

    flux_grid =[]

    for Z in models:

        file = str(path)+"bc2003_hr_xmiless_"+str(Z)+"_kroup_ssp.ised_ASCII"

        flux_models = np.genfromtxt(file, skip_header=7, skip_footer=12, dtype=float)


        flux_grid.append(flux_models[:,1:13217])

    #print(flux_grid.shape)

    return ages, wavelengths, flux_grid

#new_wavs = np.arange(1000., 70000., 10.)
#new_fluxes = spectres.spectres(new_wavs, lambda_[1:], fluxes)



#print(len(flux_grid))
#import pickle

#arr = {'ages': ages, 'wavelengths': wavelengths, 'fluxes':flux_grid}

#print(arr)
#pickle.dump( arr, open( "awf_spec.p", "wb" ))

#file = pickle.load( open( "awf_spec.p", "rb" ) )
#ages, waves,flux_grid = spectrum()
#ages = np.array(ages)
#print(file['fluxes'])
#models = flux_grid

#print(models[0:5])

#fits.writeto("awf_spec.p", arr)
