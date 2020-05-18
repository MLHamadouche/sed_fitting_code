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

mass_norm = np.loadtxt("/Users/PhDStuff/mass_normalization_bc03.txt", usecols=[10])

#filter_list =["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]
#new_filter_list = ['cfht_U.txt', 'subaru_B.txt', 'subaru_V.txt', 'subaru_R.txt', 'subaru_i.txt', 'subaru_z.txt', 'subaru_newz.txt', 'subaru_nb921.txt', 'vista_Y.txt', 'wfcam_J.txt', 'wfcam_H.txt','wfcam_K.txt']

new_filter_list = ['ECDFS_B_filter.txt', 'ECDFS_H_filter.txt',  'ECDFS_I598_filter.txt', 'ECDFS_I679_filter.txt', 'ECDFS_J_filter.txt','ECDFS_U_filter.txt', 'ECDFS_z850_filter.txt', 'ECDFS_CH1_filter.txt', 'ECDFS_I484_filter.txt',  'ECDFS_I624_filter.txt' , 'ECDFS_I738_filter.txt' ,'ECDFS_K_filter.txt' , 'ECDFS_V606_filter.txt',
'ECDFS_CH2_filter.txt' , 'ECDFS_I527_filter.txt' , 'ECDFS_I651_filter.txt',  'ECDFS_I767_filter.txt' , 'ECDFS_R_filter.txt' , 'ECDFS_Y_filter.txt']

filter_curves = pc.load_filter_files(new_filter_list)

eff_wavs = pc.calc_eff_wavs(filter_curves)

######## LOAD CATALOG DATA ##########################
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_zphot_test/ECDFS_zphot_training_phot.fits').to_pandas()
#objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
objects = np.array('ECDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))

data = []

for i in range(len(objects)):
    data.append(objects[i])

print(len(data))
ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)

####### FLUX CONVERSION #######################################
fluxes_obs = cf.conversion_func(fluxes_obs_raw, np.array(eff_wavs))
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, np.array(eff_wavs))

######## INITIALISE ARRAYS ###############################################
best_chisq = np.ones(len(data))
best_chisq*=np.inf
best_redshift = np.ones(len(data))
best_dust = np.ones(len(data))
best_ages = np.ones(len(data))
formed_mass = np.zeros(len(data))
stellar_mass = np.zeros(len(data))

time_start=time.time()
######### LOAD BC03 MODELS #############################################
file = pickle.load( open( "awf_spec.p", "rb" ) )
ages = file['ages']
flux_grid = file['fluxes']
waves = file['wavelengths']
models = flux_grid

######## SET PARAMETER RANGES #######################################
redshifts = np.arange(0.001, 6.201, 0.1)
dust_att = np.arange(0.,2.501, 0.1)
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
        for a in range(103,180,2):
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
                    age_ind = a
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
file =  "ECDFS_test_new_catalogue.fits"
hdu.writeto(file)


"""
for object in range(len(data)):
    flux_best_model = models[4][age_ind,:]

    k_lam = dusty.dust_masks(waves)
    flux_best_model *=10**(-0.4*best_dust[object]*k_lam)
    new_fluxes = pc.Photometry(waves, filter_curves, eff_wavs, best_redshift[object], flux_best_model).photometry()
    flux_best_model_plot = new_fluxes*bestest_mass[object]
    #f_lam_model = np.expand_dims(flux_best_model, axis=0)/np.expand_dims(lum_area,axis=1)
    #plt.plot(waves*(1+best_redshift[object]), f_lam_model)
    plt.scatter(eff_wavs*(1+best_redshift[object]), flux_best_model_plot, color="blue", zorder=3)
    plt.scatter(eff_wavs*(1+best_redshift[object]), fluxes_obs[object], color="red", zorder=2)
    plt.errorbar(eff_wavs*(1+best_redshift[object]), fluxes_obs[object], yerr = fluxerrs_obs[object], label='f_errors', ls=" ")
    plt.xlim(0,50000*(1+best_redshift[object]))
    plt.savefig(str(data[object])+'.png')
    plt.close()
"""
