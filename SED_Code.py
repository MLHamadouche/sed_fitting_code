import numpy as np
import matplotlib.pyplot as plt
import PhotometryClass as pc
import load_data as ld
import conversion_function as cf
import dust as dusty
import pandas as pd
import func_chisq as chisq
import pickle
import eff_wavs_filter_function as ewavs
import time
from astropy.table import Table
from astropy.io import fits


filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]

eff_wavs = ewavs.filter_wavs()
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

data = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
print(data)
data = data[0:150]

ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)

fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

best_chisq = np.ones(len(data))
best_chisq*=np.inf
best_redshift = np.ones(len(data))
best_dust = np.ones(len(data))
best_ages = np.ones(len(data))
bestest_mass = np.ones(len(data))

time_start=time.time()

file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
models = flux_grid

redshifts = np.arange(1.001, 6.201, 0.1)
dust_att = np.arange(0.1,1.901,0.1)
#103 index is 40 Million Years, 181 is 10Gyrs
total_models = ((181-103)/3)*len(redshifts)*len(dust_att)
print(f'total no. models:{total_models}')
k_lam = dusty.dust_masks(waves)

filter_curves = pc.load_filter_files(filter_list)

for z in range(len(redshifts)):
    redshift = redshifts[z]
    for d in range(len(dust_att)):
        A_v = dust_att[d]
        for a in range(103,180,3):
            model_flux = np.copy(models[4][a,:])
            time_model_start = time.time()
            no_models_done = a + d*((181-103)/3) + len(dust_att)*((181-103)/3)*z
            #if not no_models_done % 1000:
            print("Tried", no_models_done, "/", total_models, " models")

            model_flux *= 10**(-0.4*A_v*k_lam)

            #_fluxes=pc.Photometry(waves, filter_list, redshift, model_flux)
            new_fluxes = pc.Photometry(waves, filter_curves, eff_wavs, redshift, model_flux).photometry()

            best_mass = np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)),axis =1)/np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1)

            model = np.expand_dims(new_fluxes,axis =0)*np.expand_dims(best_mass, axis=1)

            diffs= model-fluxes_obs

            chisq = np.sum((diffs**2)/(fluxerrs_obs**2), axis=1)

            for m in range(len(data)):
                if chisq[m] < best_chisq[m]:
                    best_chisq[m]=chisq[m]
                    bestest_mass[m]=best_mass[m]
                    best_redshift[m]=redshift
                    best_ages[m]=ages[a]
                    best_dust[m]=A_v

            time_model_end = time.time() - time_model_start
            print(f'time model end: {time_model_end}')

time_end = time.time() - time_start
print(f'time end: {np.round(time_end/60, 3)} mins')

col1 = fits.Column(name='target', format='10A', array=data)
col2 = fits.Column(name='redshift', format='E', array=best_redshift)
col3 = fits.Column(name='age', format='E',  array=best_ages)
col4 = fits.Column(name='mass', format='E',  array=bestest_mass)
col5 = fits.Column(name='dust', format='E',  array=best_dust)
col6 = fits.Column(name='best chisq', format='E', array=best_chisq)

hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])

hdu.writeto("new_test150_catalogue.fits")
