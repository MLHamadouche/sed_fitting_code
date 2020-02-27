import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table

def read_data(data_row_index):
    table = Table.read('/Users/massissiliahamadouche/Downloads/VANDELS_CDFS_HST_PHOT_v1.0.fits')
    df = table.to_pandas()
    df.index = df['ID']

    data_row = df.index[data_row_index]


    print(data_row)
    print(data_row_index)

    flux_cols = ['CH2_flux', 'Ks_HAWKI_flux','Ks_ISAAC_flux','CH1_flux','VIMOS_U_flux','F098M_flux','F105W_flux','F125W_flux','F160W_flux','F435W_flux','F606W_flux','F775W_flux','F814W_flux','F850LP_flux']
    flux_errs_cols = ['CH2_err', 'Ks_HAWKI_err','Ks_ISAAC_err','CH1_err','VIMOS_U_err','F098M_err','F105W_err','F125W_err','F160W_err', 'F435W_err','F606W_err', 'F775W_err','F814W_err', 'F850LP_err']

    table_data = df.loc[data_row, flux_cols].values
    table_errors = df.loc[data_row, flux_errs_cols].values

    for i in range(table_data.shape[0]):
        if table_errors[i] < table_data[i]/20.0:
            table_errors[i] = table_data[i]/20.

        if (table_errors[i] < 0.) or (table_data[i] < 0.):
            table_errors[i] = 9.9*10**99
            table_data[i] = 0.

    return table_data, table_errors
#print(table_errors)
#print(np.log10(31622776601.683792))

#print(10**10.5)
