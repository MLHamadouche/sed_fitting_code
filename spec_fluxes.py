import numpy as np
import matplotlib.pyplot as plt
import spectres
import time


gen_file = "/Users/massissiliahamadouche/bc03/Miles_Atlas/Kroupa_IMF/bc2003_hr_xmiless_m62_kroup_ssp.ised_ASCII"


age_bc03 = np.loadtxt(gen_file, max_rows = 1, dtype=float)

ages = np.array(age_bc03[1:222])


lambda_bc03 = np.genfromtxt(gen_file, skip_header=6, skip_footer=233, dtype=float)
lambda_ = np.array(lambda_bc03)

wavelengths = lambda_[1:]
#new_wavs = np.arange(1000., 70000., 10.)
#new_fluxes = spectres.spectres(new_wavs, lambda_[1:], fluxes)

def spectrum():
    models = ["m22","m32", "m42","m52", "m62", "m72", "m82"]

    flux_grid =[]

    for Z in models:

        file = "/Users/massissiliahamadouche/bc03/Miles_Atlas/Kroupa_IMF/bc2003_hr_xmiless_"+str(Z)+"_kroup_ssp.ised_ASCII"

        flux_models = np.genfromtxt(file, skip_header=7, skip_footer=12, dtype=float)


        flux_grid.append(flux_models[:,1:13217])

    #print(flux_grid.shape)

    return ages, wavelengths, flux_grid

#ages, waves, flux_grid = spectrum()
#print(flux_grid[0][30][-1])
#print(len(flux_grid))
