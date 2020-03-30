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
from astropy.cosmology import FlatLambdaCDM
import load_data as ld

eff_wavs = ewavs.filter_wavs()
data = ['CDFS000005MASTER','CDFS000006SELECT','CDFS000007MASTER','CDFS000008MASTER','CDFS000009MASTER']
ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)
#ID, fluxes_obs_raw, fluxerrs_obs_raw = gd.load_data(2918)
#print(ID, fluxes_obs_raw, fluxerrs_obs_raw)

#ID, fluxes_obs_raw, fluxerrs_obs_raw = gd.read_data(2913)
#print("DATA", fluxes_obs_raw, "ERRS", fluxerrs_obs_raw)

fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

print(f'fluxes_obs.shape:{fluxes_obs.shape}')

#redshifts = np.arange(0.9, 1.16, 0.01)

redshifts = np.arange(0.2, 3.6, 0.2)
#redshifts = np.expand_dims(redshifts, axis = 1)*np.ones(5)

best_chisq = np.ones(5)
best_chisq*=np.inf
print(f'best_chisq {best_chisq}')
print(f'best_chisq.shape {best_chisq.shape}')
best_redshift = np.ones(5)
best_dust = np.ones(5)
best_ages = np.ones(5)

time_start = time.time()
import pickle

file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
models = flux_grid

print(f'models[4].shape:{models[4].shape}')

dust_att = np.arange(0.1,3.6,0.5)
#dust_att = np.expand_dims(dust_att, axis = 1)*np.ones(5)

total_models = 41*len(redshifts)*len(dust_att)
print(f'total no. models:{total_models}')

for z in range(len(redshifts)):
    redshift = redshifts[z]
    print(f'redshift start of loop: {redshift}')
    for d in range(len(dust_att)):
        A_v = dust_att[d]
        for a in range(110,150,2):
            model_flux = np.copy(models[4][a,:])
            time0 = time.time()

            no_models_done = a + d*41 + len(dust_att)*41*z
            #if not no_models_done % 1000:
            print("Tried", no_models_done, "/", total_models, " models")
            #model_flux = np.expand_dims(model_flux, axis = 0)
            k_lam = dusty.dust_masks(waves)
            #print(f'k_lam shape: {k_lam.shape}')

            model_flux *= 10**(-0.4*A_v*k_lam)
            #new_fluxes = np.expand_dims(np.zeros(5), axis =0)
            new_fluxes = pf.photometry(waves, model_flux, redshift)

            print(f'new fluxes:{new_fluxes.shape}')
            #print(new_fluxes.shape)
            print(f'sum test: {np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)), axis=1).shape}')
            print(f'sum test bottom: {np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1).shape}')
            best_mass = np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)),axis =1)/np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1)
            print(f'best mass shape {best_mass.shape}')
            print(f'best mass {best_mass}')

            model = np.expand_dims(new_fluxes,axis =0)*np.expand_dims(best_mass, axis=1)

            print(f'model shape:{model.shape}')

            #print(f'model fluxes w/ dust:{model}')
            diffs= model-fluxes_obs
            #print(f'fluxes: {fluxes_obs}')
            #print(f'fluxerrs: {fluxerrs_obs}')
            #print(f'diffs: {diffs}')
            chisq = np.sum((diffs**2)/(fluxerrs_obs**2), axis=1)
            print(f'chisq {chisq}')
            print(f'chisq_shape {chisq.shape}')
            #chisq_arr = []
            age_arr = []
            #mass_arr = []
            redshift_arr = []
            dust_arr = []
            for m in range(5):
                if chisq[m] < best_chisq[m]:
                    best_chisq=chisq
                    best_redshift[m]=redshift
                    best_ages[m]=a
                    bestest_mass=best_mass
                    best_dust[m]=A_v
                    print(redshift, a, A_v, chisq)
            age_arr.append(best_ages)
            dust_arr.append(best_dust)
            redshift_arr.append(best_redshift)
            chisq_arr=best_chisq
            #age_arr.append(best_ages)
            mass_arr=bestest_mass
            #dust_arr.append(best_dust)
            #redshift_arr.append(best_redshift)
            print(f'time phot in loop taken: {np.round(time.time() - time0, 3)}')
time_end = time.time() - time_start
print(f'time end: {np.round(time_end/60, 3)} mins')

print(f'chisq_arr: {chisq_arr}')
print(f'mass_arr: {mass_arr}')
print(f'age_arr: {age_arr}')
print(f'redshift_arr: {redshift_arr}')
print(f'dust_arr: {dust_arr}')
print(age_arr[0][0])
print(best_ages[np.int(0)])
print(f'best_redshift: {best_redshift}, bestest_mass: {bestest_mass}, best_dust: {best_dust},best_age: {best_ages}, best_chisq: {best_chisq}')


for object in range(5):

    #print(f'chisq={best_chisq[object]}, best redshift ={best_redshift[object]}, best_age={best_age[object], ages[best_age[object]]}, best_mass={bestest_mass[object], np.log10(bestest_mass[object])}, best_dust={best_dust[object]}')# best_mass={best_mass}')#best_metal={best_metal},
    #print(f'chisq={best_chisq}, best points are : z: {best_redshift}, age:{best_age}, metallicity: {best_metal}')
    flux_best_model = models[4][int(best_ages[np.int(object)]),:]

    k_lam = dusty.dust_masks(waves)
    flux_best_model *=10**(-0.4*best_dust[object]*k_lam)
    #waves, best_model = sf.spectrum(best_metal)
    flux_best_model_plot = pf.photometry(waves, flux_best_model, best_redshift[object])

    #eff_wavs1, new_fluxes, wavelengths, f_spec_model = pf.photometry(waves, best_model, best_age, best_redshift, mass)

    #plt.scatter(eff_wavs, flux_best_model_plot*bestest_mass, color="blue", zorder=3)

    plt.scatter(eff_wavs, flux_best_model_plot*bestest_mass[object], color="blue", zorder=3)

    plt.scatter(eff_wavs, fluxes_obs[object], color="red", zorder=3)
    plt.errorbar(eff_wavs, fluxes_obs[object], yerr = fluxerrs_obs[object], label='f_errors', ls=" ")
    plt.ylim(-1.5*10**-17, 1.5*10**-17)
    plt.xlim(0,50000)
    plt.savefig(str(object)+'.png')
    plt.show()


"""
    for a in range(147,161,1):
        model_flux = np.copy(models[4][a,:])


        time0 = time.time()

        no_models_done = a + d*16.9 + len(dust_att)*16.9*z
        #no_models_done = d + a*len(dust_att) + (200)*z*len(dust_att)
        #if not no_models_done % 1000:
        print("Tried", no_models_done, "/", total_models, " models")

        k_lam = dusty.dust_masks(waves)
        #print(model_flux)
        model_flux *= 10**(-0.4*A_v*k_lam)
        #print(model_flux)

        new_fluxes = pf.photometry(waves, model_flux, redshift)

        #new_fluxes *= 10**(0.4*0.28*0.44)*(A_v/A_lam)
        #print(f'new fluxes:{new_fluxes}')
        best_mass = np.sum((new_fluxes*fluxes_obs)/(fluxerrs_obs**2))/np.sum((new_fluxes**2)/(fluxerrs_obs**2))

"""
