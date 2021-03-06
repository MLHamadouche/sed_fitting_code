import numpy as np
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt

#flux_cols = ['CH2_flux', 'Ks_HAWKI_flux','Ks_ISAAC_flux','CH1_flux','VIMOS_U_flux','F098M_flux','F105W_flux','F125W_flux','F160W_flux','F435W_flux','F606W_flux','F775W_flux','F814W_flux','F850LP_flux']
#flux_errs_cols = ['CH2_err', 'Ks_HAWKI_err','Ks_ISAAC_err','CH1_err','VIMOS_U_err','F098M_err','F105W_err','F125W_err','F160W_err', 'F435W_err','F606W_err', 'F775W_err','F814W_err', 'F850LP_err']
"""
flux_cols  = ['U_flux','B_flux','V_flux','R_flux','i_flux','z_flux','newz_flux','nb921_flux','Y_flux','J_flux','H_flux','K_flux']
flux_errs_cols = ['U_err','B_err','V_err','R_err','i_err','z_err','newz_err','nb921_err','Y_err','J_err','H_err','K_err']
#catalog_file = Table.read('/massi_zphot_test/UDS_zphot_training_phot.fits').to_pandas()
"""

flux_cols = [ 'U_flux', 'B_flux', 'I484_flux' , 'I527_flux' ,'I598_flux',  'V606_flux', 'I624_flux' , 'I651_flux' ,'R', 'I679'  ,'I738' ,'I767' ,'z850_flux',  'Y_flux' , 'J_flux' , 'H_flux' , 'K_flux' , 'CH1_flux' , 'CH2_flux']
flux_errs_cols = ['U_err', 'B_err' , 'I484_err'  ,'I527_err','I598_err' ,'V606_err', 'I624_err' ,'I651_err',  'R_err' , 'I679_err', 'I738_err',  'I767_err',  'z850_err' ,  'Y_err', 'J_err' , 'H_err', 'K_err', 'CH1_err', 'CH2_err']


catalog_file = Table.read('/Users/massi_zphot_test/ecdfs_zphot_validation_phot_zspec.fits').to_pandas()


def load_catalog_data(data_array):
    catalog = pd.DataFrame(catalog_file)
    #catalog_IDs = catalog_file['ID']
    ind = catalog.set_index('ECDFS' + catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))
    #print(f'ind: {ind.index}')
    #cat_ind = ind.index
    #catalog_IDs.index = 'CDFS' + catalog_file['ID'].astype(str).str.pad(6, side="left", fillchar="0")#+ catalog_file['CAT'].str.decode("utf-8")
    #print("catalog.catalog_IDs", catalog.catalog_IDs)
    #cat = catalog_IDs.reset_index
    #print("cat[data_array]:", cat[data_array])
    #_cat = cat[data_array]
    #print("_cat:", _cat)
    #cat_data = catalog_file.loc.[_cat, flux_cols].values
    cat_data = ind.loc[data_array, flux_cols].values
    #print(cat_data)
    cat_errs = ind.loc[data_array,flux_errs_cols].values
    #print(cat_data[0].shape[0])

    for m in range(len(data_array)):
        for i in range(cat_data[m].shape[0]):
            if cat_errs[m][i] < cat_data[m][i]/20.0:
                cat_errs[m][i] = cat_data[m][i]/20.

            if (cat_errs[m][i] < 0.) or (cat_data[m][i] < 0.):
                cat_errs[m][i] = 9.9*10**99
                cat_data[m][i] = 0.

    return ind.index, cat_data, cat_errs
