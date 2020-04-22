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
import bagpipes as pipes
#from tqdm import tqdm

filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]

eff_wavs = ewavs.filter_wavs()
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))

#print(data)
data = []

for i in range(len(objects)):
    data.append(objects[i])
#print(data)


ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)

fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)

input()
best_chisq = np.ones(len(data))
best_chisq*=np.inf
best_redshift = np.ones(len(data))
best_dust = np.ones(len(data))
best_ages = np.ones(len(data))
bestest_mass = np.zeros(len(data))

time_start=time.time()

file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
models = flux_grid

redshifts = np.arange(1.001, 6.201, 0.05)
dust_att = np.arange(0.,2.501,0.1)
#103 index is 40 Million Years, 181 is 10Gyrs
total_models = ((181-103)/2)*len(redshifts)*len(dust_att)
print(f'total no. models:{total_models}')
k_lam = dusty.dust_masks(waves)

filter_curves = pc.load_filter_files(filter_list)

igm = pipes.models.igm(waves)


for z in range(len(redshifts)):
    redshift = redshifts[z]
    for d in range(len(dust_att)):
        A_v = dust_att[d]
        for a in range(103,180,2):
            model_flux = np.copy(models[4][a,:])
            time_model_start = time.time()
            no_models_done = a + d*((181-103)/2) + len(dust_att)*((181-103)/2)*z
            #if not no_models_done % 1000:
            print("Tried", no_models_done, "/", total_models, " models")

            model_flux *= 10**(-0.4*A_v*k_lam)
            igm_trans = igm.trans(redshift)
            model_flux*=igm_trans
            #print(f'igm_trans:{igm_trans}')
            #_fluxes=pc.Photometry(waves, filter_list, redshift, model_flux)
            new_fluxes = pc.Photometry(waves, filter_curves, eff_wavs, redshift, model_flux).photometry()
            #print(new_fluxes)

            #print(new_fluxes)
            best_mass = np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)),axis =1)/np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1)
            #print(f'best_mass: {best_mass}')

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

            #time_model_end = time.time() - time_model_start
            #print(f'time model end: {time_model_end}')



time_end = time.time() - time_start
print(f'time end: {np.round(time_end/60, 3)} mins')
RA = ross_objects['RA']
DEC = ross_objects['DEC']

col1 = fits.Column(name='target', format='10A', array=data)
col2 = fits.Column(name='redshift', format='E', array=best_redshift)
col3 = fits.Column(name='age', format='E',  array=best_ages)
col4 = fits.Column(name='mass', format='E',  array=bestest_mass)
col5 = fits.Column(name='dust', format='E',  array=best_dust)
col6 = fits.Column(name='best chisq', format='E', array=best_chisq)
col7 = fits.Column(name='RA', format='E', array=RA)
col8 = fits.Column(name='DEC', format='E', array=DEC)

hdu = fits.BinTableHDU.from_columns([col1, col7, col8, col2, col3, col4, col5, col6])
file =  "masstest2_ra_dec_catalogue.fits"
hdu.writeto(file)
