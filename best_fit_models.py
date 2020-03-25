import numpy as np
import matplotlib.pyplot as plt
import conversion_function as cf
import time
import phot_filts as pf
import eff_wavs_filter_function as ewavs
import dust as dusty
from astropy.io import fits
from astropy.table import Table
import load_data as ld
import pickle

eff_wavs = ewavs.filter_wavs()

# READING IN DATA FROM CATALOG FILE ############################################
data = np.arange(5)
ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)

fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

print(f'fluxes_obs.shape:{fluxes_obs.shape}')

# READING IN BC03 MODELS #######################################################
file = pickle.load( open( "awf_spec.p", "rb" ) )
age_ind = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']

model_test = np.ones(5)

models = flux_grid


print(f'models[4].shape:{models[4].shape}')

dust_att = np.arange(0.1,4.5,0.5)
#print(dust_att)

z_ = np.arange(0.2, 3.6, 0.2)
#print(z_)

best_chisq = np.ones(5)
best_chisq*=np.inf
#print(f'best_chisq {best_chisq}')
#print(f'best_chisq.shape {best_chisq.shape}')


#ages = np.arange(1,222, 1)
#print(ages)
model_grid = np.expand_dims(np.ones_like(age_ind), axis =0)*age_ind
#print(f'one: {model_grid.shape}')
model_grid =  np.expand_dims(z_, axis = 1)*model_grid
model_grid  = np.expand_dims(np.ones_like(model_grid), axis = 0) * np.expand_dims(np.ones_like(dust_att), axis = [1,2])
#model_grid = np.expand_dims(model_grid, axis = 0)*dust_att

#print(f'three: {model_grid.shape}')

grid = np.copy(model_grid)

dust = np.ones(5)
#print(dust)
dust = np.expand_dims(dust, axis = 0)*np.expand_dims(dust_att, axis = 1)
redshifts = np.ones(5)
redshifts = np.expand_dims(redshifts, axis = 0)* np.expand_dims(z_, axis = 1)
ages= np.ones(5)
ages = np.expand_dims(ages, axis = 0)*np.expand_dims(age_ind, axis =1)

print(f'redshift shape before loop: {redshifts.shape}')

print(f'dust before loop: {dust.shape}')

for z in range(len(z_)):
    redshift = redshifts[:,z]
    for d in range(len(dust_att)):
        A_v = dust[:,d]
        for a in range(len((100,180))):

            model_flux = np.copy(models[4][a,:])

            model_flux_Z  = np.expand_dims(model_flux, axis =0)

            print(f'modelflux shape: {model_flux_Z.shape}')
            #print(f'dust shape: {A_v.shape}')

            time0 = time.time()

            #no_models_done = a + d*221 + 9*221*z
            #if not no_models_done % 1000:
            #print("Tried", no_models_done, "/", total_models, " models")

            k_lam = dusty.dust_masks(waves)
            print(f'A_v.shape: {A_v.shape}')
            #print(model_flux)
            model_flux_ = model_flux_Z * (10**(-0.4*np.expand_dims(A_v, axis = 1)*k_lam))
            print(f'model_flux_Z.shape: {model_flux_Z.shape}')
            #print(f'A_v.shape: {A_v.shape}')
            #print(f'redshift.shape: {redshift.shape}')
            print(f'modelflux shape: {model_flux_.shape}')
            print(f'redshift: {redshift.shape}')
            #waves = np.expand_dims(waves, axis = 0)*np.expand_dims(np.ones(5), axis = 1)
            #print(f'waves.shape: {waves.shape}')

            new_fluxes = pf.photometry(waves, model_flux_, redshift)

            print(f'waves.shape: {waves.shape}')
            print(f'new fluxes shape:{new_fluxes.shape}')
            #print(new_fluxes.shape)
            #print(fluxes_obs.shape)
            #print(fluxerrs_obs.shape)
            #print(f'sum test: {np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)), axis=1).shape}')
            #print(f'sum test bottom: {np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1).shape}')
            best_mass = np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)),axis =1)/np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1)
            print(f'best mass shape {best_mass.shape}')
            print(f'best mass {best_mass}')

            model = new_fluxes*np.expand_dims(best_mass, axis=1)
            #model = np.expand_dims(new_fluxes, axis=0)*np.expand_dims(best_mass, axis=1)

            print(f'model shape:{model.shape}')

            #print(f'model fluxes w/ dust:{model}')
            diffs= model-fluxes_obs
            #print(f'fluxes: {fluxes_obs}')
            #print(f'fluxerrs: {fluxerrs_obs}')
            #print(f'diffs: {diffs}')
            chisq = np.sum((diffs**2)/(fluxerrs_obs**2), axis=1)
            print(f'chisq {chisq}')
            print(f'chisq_shape {chisq.shape}')

            if (chisq < best_chisq).all :
                best_chisq = chisq
                best_redshift = redshift
                best_age = a
                bestest_mass = best_mass
                best_dust = A_v
                print(redshift, ages[a], A_v, chisq)

            print(f'time phot in loop taken: {np.round(time.time() - time0, 3)}')

time_end = time.time() - time_start
print(f'time end: {np.round(time_end/60, 3)} mins')
print(f'chisq={best_chisq}, best redshift ={best_redshift}, best_age={best_age, ages[best_age]}, best_mass={bestest_mass, np.log10(bestest_mass)}, best_dust={best_dust}')# best_mass={best_mass}')#best_metal={best_metal},
#print(f'chisq={best_chisq}, best points are : z: {best_redshift}, age:{best_age}, metallicity: {best_metal}')
flux_best_model = models[4][best_age, :]

k_lam = dusty.dust_masks(waves)
flux_best_model *=10**(-0.4*best_dust*k_lam)
#waves, best_model = sf.spectrum(best_metal)
flux_best_model_plot = pf.photometry(waves, flux_best_model, best_redshift)

#eff_wavs1, new_fluxes, wavelengths, f_spec_model = pf.photometry(waves, best_model, best_age, best_redshift, mass)

#plt.scatter(eff_wavs, flux_best_model_plot*bestest_mass, color="blue", zorder=3)

plt.scatter(eff_wavs, flux_best_model_plot*bestest_mass, color="blue", zorder=3)

plt.scatter(eff_wavs, fluxes_obs, color="red", zorder=3)
plt.errorbar(eff_wavs, fluxes_obs, yerr = fluxerrs_obs, label='f_errors', ls=" ")
plt.ylim(-1.5*10**-17, 1.5*10**-17)
plt.xlim(0,50000)
plt.show()
