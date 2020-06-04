import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import os
import re
#flux_cols  = ['U_flux','B_flux','V_flux','R_flux','i_flux','z_flux','newz_flux','nb921_flux','Y_flux','J_flux','H_flux','K_flux']
#flux_errs_cols = ['U_err','B_err','V_err','R_err','i_err','z_err','newz_err','nb921_err','Y_err','J_err','H_err','K_err']

def find_file(ID, extension):
    #ID = ID.split()[1]
    new_ID = re.search('\d+', ID).group()

    for root, dirs, files in os.walk('/Users/PhDStuff/sed_fitting_code'):
        if ID.startswith("CDFS"):
            files = [filename for filename in files if filename.startswith("VANDELS_CDFS") and filename.endswith(".{ext}".format(ext=extension))]

        elif ID.startswith("UDS"):
            files = [filename for filename in files if filename.startswith("VANDELS_UDS") and filename.endswith(".{ext}".format(ext=extension))]

        for filename in files:
            path = os.path.join(root, filename)
            print(path)
            data_file = Table.read(path).to_pandas()
            #print(data_file)
            catalog = pd.DataFrame(data_file)
            ID_list = catalog['ID']
            print(ID_list)
            if new_ID in ID_list:
                both = os.path.split(path)
                print(both)
                print(path)
                prefix_for_load = both[1]
                print(prefix_for_load)
                prefix = prefix_for_load.split('VANDELS_')[1]
                pre = prefix.split('_PHOT')[0]
                #print(pre)
                return path, pre
#data_array = 'CDFS247586'
path, prefix = find_file('CDFS416284', 'fits')

print(f'path: {path}, prefix:{prefix}')

cat_file = Table.read(path).to_pandas()
catalog = pd.DataFrame(cat_file)

def load_vandels(data_array):

    ind = catalog.set_index(str(prefix) + catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0")) + catalog['CAT'].astype('udf8')
    # Extract the object we want from the catalogue.
    fluxes = ind.loc[data_array, flux_cols].values

    fluxerrs = ind.loc[data_array,flux_errs_cols].values

    photometry = np.c_[fluxes, fluxerrs]

    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]

    for i in range(len(photometry)):
        if i < 10:
            max_snr = 20.

        else:
            max_snr = 10.

        if photometry[i, 0]/photometry[i, 1] > max_snr:
            photometry[i, 1] = photometry[i, 0]/max_snr

    return photometry
