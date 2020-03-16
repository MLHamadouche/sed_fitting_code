import numpy as np
import matplotlib.pyplot as plt
import conversion_function as cf
import get_data as gd
import time
import spec_fluxes as sf
import phot_filts as pf
import eff_wavs_filter_function as ewavs
import dust as dusty
from astropy.io import fits
from astropy.table import Table
import pickle


file = pickle.load( open( "awf_spec.p", "rb" ) )
#ages, waves,flux_grid = spectrum()
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']

def catalog(ID_list):


        eff_wavs = ewavs.filter_wavs()

        models = flux_grid

        redshifts = np.arange(0.5,3.5,0.1)


        dust_att = np.arange(0.5,2.0,0.1)

        total_models = 35*len(redshifts)*len(dust_att)
        print(f'total no. models:{total_models}')


        obj_ID, fluxes_obs_raw, fluxerrs_obs_raw = gd.read_data(ID_list)

        fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
        fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

        best_chisq = np.inf

        time_start = time.time()


        for z in range(len(redshifts)):
            redshift = redshifts[z]
            for d in range(len(dust_att)):
                A_v = dust_att[d]
                for a in range(100,180,1):
                    model_flux = np.copy(models[4][a,:])

                    time0 = time.time()

                    no_models_done = a + d*35 + 35*len(dust_att)*z
                    #no_models_done = d + a*len(dust_att) + (200)*z*len(dust_att)
                    #if not no_models_done % 1000:
                    print("Tried", no_models_done, "/", total_models, " models")
                    k_lam = dusty.dust_masks(waves)
                    #print(model_flux)
                    model_flux *= 10**(-0.4*A_v*k_lam)

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
                        best_dust= A_v
                        print(redshift, ages[a], chisq)



        time_end = time.time() - time_start
        print(f'time end: {np.round(time_end, 3)}')
        print(f'chisq={best_chisq}, best redshift ={best_redshift}, best_age={best_age, ages[best_age]}, best_mass={bestest_mass, np.log10(bestest_mass)}, best_dust={best_dust}')
        #print(f'chisq={best_chisq}, best points are : z: {best_redshift}, age:{best_age}, metallicity: {best_metal}')
        flux_best_model = models[4][best_age, :]
        k_lam = dusty.dust_masks(waves)
        flux_best_model *=10**(-0.4*best_dust*k_lam)
        #waves, best_model = sf.spectrum(best_metal)
        flux_best_model_plot = pf.photometry(waves, flux_best_model, best_redshift)
        #waves, best_model = sf.spectrum(best_metal)
        #flux_model = pf.photometry(waves, models[4][best_age, :], best_redshift)

        flux_best_model_plot*=bestest_mass


        #plt.scatter(eff_wavs, flux_best_model_plot, color="blue", zorder=3)
        #plt.scatter(eff_wavs, fluxes_obs, color="red", zorder=3)
        #plt.errorbar(eff_wavs, fluxes_obs, yerr = fluxerrs_obs, label='f_errors', ls=" ")
        #plt.ylim(-3*10**-18, 3*10**-18)
        #plt.xlim(0,50000)
        #plt.show()



    return obj_ID, best_redshift, bestest_mass, np.log10(bestest_mass), ages[best_age], best_dust

ID_list = [6,31]
print(ID_list)

obj_ID, best_z, best_m, log_m, best_a, best_dust = catalog(ID_list)
#input()
