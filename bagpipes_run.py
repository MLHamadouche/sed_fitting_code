import bagpipes as pipes
import numpy as np
import PhotometryClass as pc
import conversion_function as cf
from astropy.io import fits
from astropy.table import Table
import pandas as pd

burst = {}
burst["age"] = (0., 13.)           # Vary the age between 0 and 13 Gyr
burst["metallicity"] = 1.0  # Vary stellar metallicity between 0 and 5 times Solar
burst["massformed"] = (0., 13.)    # Vary the Log_10 of total mass formed from 0 to 13.

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0., 2.5)

fit_instructions = {}
fit_instructions["redshift"] = (0., 10.)
fit_instructions["burst"] = burst
fit_instructions["dust"] = dust

"""
filters = ["/Users/PhDStuff/sed_fitting_code/filters/CH2",
            "/Users/PhDStuff/sed_fitting_code/filters/HAWKI_K",
            "/Users/PhDStuff/sed_fitting_code/filters/ISAAC_Ks",
            "/Users/PhDStuff/sed_fitting_code/filters/CH1",
            "/Users/PhDStuff/sed_fitting_code/filters/VIMOS_U","/Users/PhDStuff/sed_fitting_code/filters/f098m",
            "/Users/PhDStuff/sed_fitting_code/filters/f105w","/Users/PhDStuff/sed_fitting_code/filters/f125w",
            "/Users/PhDStuff/sed_fitting_code/filters/f160w","/Users/PhDStuff/sed_fitting_code/filters/f435w","/Users/PhDStuff/sed_fitting_code/filters/f606w",
            "/Users/PhDStuff/sed_fitting_code/filters/f775w","/Users/PhDStuff/sed_fitting_code/filters/f814w", "/Users/PhDStuff/sed_fitting_code/filters/f850lp"]


#filters_list = np.loadtxt("//Users/massissiliahamadouche/anaconda3/lib/python3.7/site-packages/bagpipes/filters/massi.filt_list", dtype='str')
exp = {}                          # Tau model star formation history component
exp["age"] = 3.                   # Gyr
exp["tau"] = 0.75                 # Gyr
exp["massformed"] = 9.            # log_10(M*/M_solar)
exp["metallicity"] = 0.5
dust = {}                         # Dust component
dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
dust["Av"] = 0.2
model_components = {}                   # The model components dictionary
model_components["redshift"] = 1.0      # Observed redshift
model_components["exponential"] = exp
model_components["dust"] = dust

model = pipes.model_galaxy(model_components, filt_list = filters)

#print(filters_list)

#model = pipes.model_galaxy(model_components, filt_list=filters_list)


fig = model.plot()
fig = model.sfh.plot()

"""
#model = pipes.model_galaxy(model_components, spec_wavs=np.arange(5000., 10000., 5.))

#fig = model.plot()


#new_filter_list = [ 'ECDFS_U_filter.txt', 'ECDFS_B_filter.txt', 'ECDFS_I484_filter.txt' , 'ECDFS_I527_filter.txt' ,'ECDFS_I598_filter.txt',  'ECDFS_V606_filter.txt', 'ECDFS_I624_filter.txt' , 'ECDFS_I651_filter.txt' ,'ECDFS_R_filter.txt', 'ECDFS_I679_filter.txt'  ,'ECDFS_I738_filter.txt' ,'ECDFS_I767_filter.txt' ,'ECDFS_z850_filter.txt',  'ECDFS_Y_filter.txt' , 'ECDFS_J_filter.txt' , 'ECDFS_H_filter.txt' , 'ECDFS_K_filter.txt' , 'ECDFS_CH1_filter.txt' , 'ECDFS_CH2_filter.txt']

#catalog = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()
new_filter_list = [ '/Users/massi_zphot_test/ECDFS_filters/ECDFS_U_filter.txt', '/Users/massi_zphot_test/ECDFS_filters/ECDFS_B_filter.txt', '/Users/massi_zphot_test/ECDFS_filters/ECDFS_I484_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_I527_filter.txt' ,'/Users/massi_zphot_test/ECDFS_filters/ECDFS_I598_filter.txt',  '/Users/massi_zphot_test/ECDFS_filters/ECDFS_V606_filter.txt', '/Users/massi_zphot_test/ECDFS_filters/ECDFS_I624_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_I651_filter.txt' ,'/Users/massi_zphot_test/ECDFS_filters/ECDFS_R_filter.txt', '/Users/massi_zphot_test/ECDFS_filters/ECDFS_I679_filter.txt'  ,'/Users/massi_zphot_test/ECDFS_filters/ECDFS_I738_filter.txt' ,'/Users/massi_zphot_test/ECDFS_filters/ECDFS_I767_filter.txt' ,'/Users/massi_zphot_test/ECDFS_filters/ECDFS_z850_filter.txt',  '/Users/massi_zphot_test/ECDFS_filters/ECDFS_Y_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_J_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_H_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_K_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_CH1_filter.txt' , '/Users/massi_zphot_test/ECDFS_filters/ECDFS_CH2_filter.txt']


#filter_curves = pc.load_filter_files('/Users/massi_zphot_test/ECDFS_filters/', new_filter_list )
#print(filter_curves)
#eff_wavs = pc.calc_eff_wavs(filter_curves)


#flux_cols = ['CH2_flux', 'Ks_HAWKI_flux','Ks_ISAAC_flux','CH1_flux','VIMOS_U_flux','F098M_flux','F105W_flux','F125W_flux','F160W_flux','F435W_flux','F606W_flux','F775W_flux','F814W_flux','F850LP_flux']
#flux_errs_cols = ['CH2_err', 'Ks_HAWKI_err','Ks_ISAAC_err','CH1_err','VIMOS_U_err','F098M_err','F105W_err','F125W_err','F160W_err', 'F435W_err','F606W_err', 'F775W_err','F814W_err', 'F850LP_err']

catalog = Table.read('/Users/massi_zphot_test/ecdfs_zphot_validation_phot_zspec.fits').to_pandas()
objects = np.array('ECDFS'+ catalog['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
flux_cols = [ 'U_flux', 'B_flux', 'I484_flux' , 'I527_flux' ,'I598_flux',  'V606_flux', 'I624_flux' , 'I651_flux' ,'R', 'I679'  ,'I738' ,'I767' ,'z850_flux',  'Y_flux' , 'J_flux' , 'H_flux' , 'K_flux' , 'CH1_flux' , 'CH2_flux']
flux_errs_cols = ['U_err', 'B_err' , 'I484_err'  ,'I527_err','I598_err' ,'V606_err', 'I624_err' ,'I651_err',  'R_err' , 'I679_err', 'I738_err',  'I767_err',  'z850_err' ,  'Y_err', 'J_err' , 'H_err', 'K_err', 'CH1_err', 'CH2_err']

data = []

for i in range(len(objects)):
    data.append(objects[i])


def load_data(data_array):

    ind = catalog.set_index('ECDFS' + catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))

    fluxes = ind.loc[data_array, flux_cols].values

    fluxerrs = ind.loc[data_array,flux_errs_cols].values
    print(fluxes)

    photometry = np.c_[fluxes, fluxerrs]

    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]

    # Enforce a maximum SNR of 20, or 10 in the IRAC channels.
    for i in range(len(photometry)):
        if i < 10:
            max_snr = 20.

        else:
            max_snr = 10.

        if photometry[i, 0]/photometry[i, 1] > max_snr:
            photometry[i, 1] = photometry[i, 0]/max_snr

    return photometry


fit_cat = pipes.fit_catalogue(data, fit_instructions, load_data, spectrum_exists=False,
                              cat_filt_list=new_filter_list, make_plots=True, run="ecdfs_cat")

fit_cat.fit(verbose=True, mpi_serial=True)
