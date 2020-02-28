import numpy as np
import matplotlib.pyplot as plt
import conversion_function as cf
import get_data as gd
import time
import spec_fluxes as sf
import phot_filts as pf
import eff_wavs_filter_function as ewavs
import dust_attenuation as dusty
from astropy.io import fits
from astropy.table import Table



def catalog(ID_list):

    for ID in ID_list:

        eff_wavs = ewavs.filter_wavs()
        ages, waves,flux_grid = sf.spectrum()

        models = flux_grid

        redshifts = np.arange(0.8, 4, 0.1)


        age = np.arange(120,180,1)

        total_models = len(age)*len(redshifts)
        print(f'total no. models:{total_models}')

        obj_ID, fluxes_obs_raw, fluxerrs_obs_raw = gd.read_data(ID)

        fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
        fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

        best_chisq = np.inf

        time_start = time.time()


        for z in range(len(redshifts)):
            redshift = redshifts[z]
            for a in range(len(age)):
                model_flux = np.copy(models[4][a,:])

                time0 = time.time()

                no_models_done = a + z*len(age)
                #no_models_done = d + a*len(dust_att) + (200)*z*len(dust_att)
                #if not no_models_done % 1000:
                print("Tried", no_models_done, "/", total_models, " models")

                new_fluxes = np.copy(pf.photometry(waves, model_flux, redshift))

                #new_fluxes *= 10**(0.4*0.28*0.44)*(A_v/A_lam)
                #print(f'new fluxes:{new_fluxes}')
                best_mass = np.sum((new_fluxes*fluxes_obs)/(fluxerrs_obs**2))/np.sum((new_fluxes**2)/(fluxerrs_obs**2))

                model = new_fluxes*best_mass
                #print(f'model values:{model}')

                diffs= model-fluxes_obs
                #print(f'fluxes: {fluxes_obs}')
                #print(f'fluxerrs: {fluxerrs_obs}')
                #print(f'diffs: {diffs}')
                chisq = np.sum((diffs**2)/(fluxerrs_obs**2))
                #input()
                #print(f'chisq: {chisq}')
                #print(f'time phot in loop taken: {np.round(time.time() - time0, 3)}')
                #input()
                if chisq<best_chisq:
                    best_chisq = chisq
                    best_redshift = redshift
                    best_age = a
                    bestest_mass = best_mass
                    #best_dust_ratio = A_v/A_lam
                    print(redshift, ages[a], chisq)



        time_end = time.time() - time_start
        print(f'time end: {np.round(time_end, 3)}')
        print(f'chisq={best_chisq}, best redshift ={best_redshift}, best_age={best_age, ages[best_age]}, best_mass={bestest_mass, np.log10(bestest_mass)}')
        #print(f'chisq={best_chisq}, best points are : z: {best_redshift}, age:{best_age}, metallicity: {best_metal}')

        #waves, best_model = sf.spectrum(best_metal)
        flux_model = pf.photometry(waves, models[4][best_age, :], best_redshift)

        flux_best_model=flux_model*bestest_mass


        plt.scatter(eff_wavs, flux_best_model, color="blue", zorder=3)
        plt.scatter(eff_wavs, fluxes_obs, color="red", zorder=3)
        plt.errorbar(eff_wavs, fluxes_obs, yerr = fluxerrs_obs, label='f_errors', ls=" ")
        plt.ylim(-3*10**-18, 3*10**-18)
        #plt.xlim(0,50000)
        plt.show()



    return obj_ID, best_redshift, bestest_mass, np.log10(bestest_mass), ages[best_age]

ID_list = [6,7,8,9]
print(ID_list)

obj_ID, best_z, best_m, log_m, best_a = catalog(ID_list)
#input()
catalog_array = []
for ID in ID_list:
    catalog_array.append([ID, best_z[ID], best_m[ID], log_m[ID], best_a[ID]])
print(catalog_array)
