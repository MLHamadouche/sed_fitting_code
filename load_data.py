import numpy as np
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt

#data_array = np.arange(10)
flux_cols = ['CH2_flux', 'Ks_HAWKI_flux','Ks_ISAAC_flux','CH1_flux','VIMOS_U_flux','F098M_flux','F105W_flux','F125W_flux','F160W_flux','F435W_flux','F606W_flux','F775W_flux','F814W_flux','F850LP_flux']
flux_errs_cols = ['CH2_err', 'Ks_HAWKI_err','Ks_ISAAC_err','CH1_err','VIMOS_U_err','F098M_err','F105W_err','F125W_err','F160W_err', 'F435W_err','F606W_err', 'F775W_err','F814W_err', 'F850LP_err']
catalog_file = Table.read('/Users/massissiliahamadouche/Downloads/VANDELS_CDFS_HST_PHOT_v1.0.fits').to_pandas()

def load_catalog_data(data_array):

    IDs = (str(i) for i in data_array)

    catalog_file.index = catalog_file['ID']
    _cat = catalog_file.index

    cat_ind = catalog_file.loc[_cat[data_array], flux_cols]
    cat_err_ind = catalog_file.loc[_cat[data_array], flux_errs_cols]

    return _cat[data_array], cat_ind, cat_err_ind

data = np.arange(5)

ID, flux, errs = load_catalog_data(data)

print(f' ID:{ID}, flux: {flux}, errs: {errs}')
