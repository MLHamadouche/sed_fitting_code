import numpy as np
import matplotlib.pyplot as plt
import PhotometryClass as pc
import load_data as ld
import conversion_function as cf
import dust as dusty
import pandas as pd
#import spec_fluxes as sf
import time
from astropy.table import Table
from astropy.io import fits
import bagpipes as pipes
import pickle
import os

if not os.path.exists("awf_spec.p"):
    ages, waves, flux_grid = sf.spectrum()
    arr = {'ages': ages, 'wavelengths': waves, 'fluxes': flux_grid}
    #print(arr)
    pickle.dump( arr, open( "awf_spec.p", "wb" ))

mass_norm = np.loadtxt("mass_normalization_bc03.txt", usecols=[10])
#filters must be in ascending wavelength order

#filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]
#print(len(new_filter_list))
#put in wavelength order
#new_filter_list = ['cfht_U.txt', 'subaru_B.txt', 'subaru_V.txt', 'subaru_R.txt', 'subaru_i.txt', 'subaru_z.txt', 'subaru_newz.txt', 'subaru_nb921.txt', 'vista_Y.txt', 'wfcam_J.txt', 'wfcam_H.txt','wfcam_K.txt']
new_filter_list = [ 'ECDFS_U_filter.txt', 'ECDFS_B_filter.txt', 'ECDFS_I484_filter.txt' , 'ECDFS_I527_filter.txt' ,'ECDFS_I598_filter.txt',  'ECDFS_V606_filter.txt', 'ECDFS_I624_filter.txt' , 'ECDFS_I651_filter.txt' ,'ECDFS_R_filter.txt', 'ECDFS_I679_filter.txt'  ,'ECDFS_I738_filter.txt' ,'ECDFS_I767_filter.txt' ,'ECDFS_z850_filter.txt',  'ECDFS_Y_filter.txt' , 'ECDFS_J_filter.txt' , 'ECDFS_H_filter.txt' , 'ECDFS_K_filter.txt' , 'ECDFS_CH1_filter.txt' , 'ECDFS_CH2_filter.txt']

#dir_path_filters = os.path.realpath('massi_zphot_test/ECDFS_filters/')

filter_curves = pc.load_filter_files('massi_zphot_test/ECDFS_filters/', new_filter_list )

eff_wavs = pc.calc_eff_wavs(filter_curves)

######## LOAD CATALOG DATA ##########################

ross_objects = Table.read('massi_zphot_test/ECDFS_zphot_training_phot.fits').to_pandas()
#objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
objects = np.array('ECDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))

data = []

for i in range(len(objects)):
    data.append(objects[i])

#print(len(data))
ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)
#print(len(fluxes_obs_raw[1]))

#func that takes those above and plots the data for that object

####### FLUX CONVERSION #######################################
fluxes_obs = cf.conversion_func(fluxes_obs_raw, np.array(eff_wavs))
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, np.array(eff_wavs))
#print(fluxes_obs[1])

plt.scatter(eff_wavs, fluxes_obs[1206], color ='mediumorchid', marker='o', linestyle='None')
plt.errorbar(eff_wavs, fluxes_obs[1206], yerr=fluxerrs_obs[1206], linestyle='None', color='k')
ax = plt.gca()
# recompute the ax.dataLim
ax.relim()
# update ax.viewLim using the new dataLim
ax.autoscale_view()
plt.draw()
plt.savefig('phot_plots/test_phot.png')
plt.close()

input()
######## INITIALISE ARRAYS ###############################################
best_chisq = np.ones(len(data))
best_chisq*=np.inf
best_redshift = np.ones(len(data))
best_dust = np.ones(len(data))
best_ages = np.ones(len(data))
formed_mass = np.zeros(len(data))
stellar_mass = np.zeros(len(data))
age_ind = np.zeros(len(data))

time_start=time.time()
######### LOAD BC03 MODELS #############################################
file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
models = flux_grid

######## SET PARAMETER RANGES #######################################
redshifts = np.arange(0.001, 5.201, 0.05)
dust_att = np.arange(0., 2.5, 0.1)
#103 index is 40 Million Years, 181 is 10Gyrs
total_models = ((181-103)/2)*len(redshifts)*len(dust_att)
print(f'total no. models:{total_models}')

######### DUST AND IGM MODEL #######################################
k_lam = dusty.dust_masks(waves)
igm = pipes.models.igm(waves)

######### MODELS LOOP ############################################

for z in range(len(redshifts)):
    redshift = redshifts[z]
    for d in range(len(dust_att)):
        A_v = dust_att[d]
        for a in range(103,181,2):
            model_flux = np.copy(models[4][a,:])
            time_model_start = time.time()
            no_models_done = a + d*((181-103)/2) + len(dust_att)*((181-103)/2)*z
            #if not no_models_done % 1000:
            print("Tried", no_models_done, "/", total_models, " models")

            model_flux *= 10**(-0.4*A_v*k_lam)
            igm_trans = igm.trans(redshift)
            model_flux*=igm_trans

            new_fluxes = pc.Photometry(waves, filter_curves, eff_wavs, redshift, model_flux).photometry()

            best_mass=np.sum(((new_fluxes*fluxes_obs)/(fluxerrs_obs**2)),axis =1)/np.sum((new_fluxes**2)/(fluxerrs_obs**2), axis=1)
            #print(f'best_mass: {best_mass}')
            model = np.expand_dims(new_fluxes,axis =0)*np.expand_dims(best_mass, axis=1)

            diffs= model-fluxes_obs

            chisq = np.sum((diffs**2)/(fluxerrs_obs**2), axis=1)

            for m in range(len(data)):
                if chisq[m] < best_chisq[m]:
                    best_chisq[m]=chisq[m]
                    best_ages[m]=ages[a]
                    age_ind[m] = a
                    formed_mass[m]=best_mass[m]
                    stellar_mass[m]=best_mass[m]*mass_norm[a]
                    best_redshift[m]=redshift
                    best_dust[m]=A_v


#pipes.models.making.igm_inoue2014.test()
time_end = time.time() - time_start
print(f'time end: {np.round(time_end/60, 3)} mins')

####### SAVE DATA AS FITS CATALOGUE ##########################################

#RA = ross_objects['RA']
#DEC = ross_objects['DEC']
col1 = fits.Column(name='target', format='10A', array=data)
col2 = fits.Column(name='redshift', format='E', array=best_redshift)
col3 = fits.Column(name='age', format='E',  array=best_ages)
col4 = fits.Column(name='formed_mass', format='E', array=formed_mass)
col9 = fits.Column(name='stellar_mass', format='E', array=stellar_mass)
col5 = fits.Column(name='dust', format='E',  array=best_dust)
col6 = fits.Column(name='best chisq', format='E', array=best_chisq)
#col7 = fits.Column(name='RA', format='E', array=RA)
#col8 = fits.Column(name='DEC', format='E', array=DEC)

#hdu = fits.BinTableHDU.from_columns([col1, col7, col8, col2, col3, col4, col9, col5, col6 ])
hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col9, col5, col6 ])
file =  "ECDFS_plots_new_catalogue.fits"
hdu.writeto(file)

matplotlib.rcParams['font.family'] = "AppleMyungjo"
#uncomment to plot the model photometry and  for each object
"""
for object in range(len(data)):
    age_index = np.int(age_ind[object])
    #print(age_index)
    new_model_flux = np.copy(models[4][age_index])

    new_model_flux *= 10**(-0.4*best_dust[object]*k_lam)
    igm_trans = igm.trans(best_redshift[object])
    new_model_flux*=igm_trans


    new_fluxes = pc.Photometry(waves, filter_curves, eff_wavs, best_redshift[object], new_model_flux).photometry()

    new_fluxes*=stellar_mass[object]

    plt.scatter(eff_wavs, fluxes_obs[object], color ='mediumorchid', marker='o', linestyle='None')
    plt.scatter(eff_wavs, new_fluxes, color ='steelblue', marker='o', linestyle='None')
    plt.errorbar(eff_wavs, fluxes_obs[object], yerr=fluxerrs_obs[object], linestyle='None', color='k')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'F$_{\lambda}$ ($erg s^{-1} cm^{-2} \AA^{-1}$)')
    ax = plt.gca()
    # recompute the ax.dataLim
    ax.relim()
    # update ax.viewLim using the new dataLim
    ax.autoscale_view()
    plt.draw()
    plt.savefig('phot_plots/'+str(data[object])+'_phot.png')
    plt.close()
"""
