import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import os
import re

def find_file(ID, extension):
    new_ID = re.search('\d+', ID).group()
    new_ID = new_ID.lstrip('0')
    print(new_ID)
    for root, dirs, files in os.walk('/Users/PhDStuff/sed_fitting_code'):
        if ID.startswith("CDFS"):
            files = [filename for filename in files if filename.startswith("VANDELS_CDFS") and filename.endswith(".{ext}".format(ext=extension))]

        elif ID.startswith("UDS"):
            files = [filename for filename in files if filename.startswith("VANDELS_UDS") and filename.endswith(".{ext}".format(ext=extension))]

        for filename in files:
            path = os.path.join(root, filename)

            data_file = Table.read(path).to_pandas()
            catalog = pd.DataFrame(data_file)
            ID_list = catalog['ID'].astype(str)

            if new_ID in ID_list.values:
                print(new_ID)
                both = os.path.split(path)
                print(both)
                print(path)
                prefix_for_load = both[1]
                print(prefix_for_load)
                prefix = prefix_for_load.split('VANDELS_')[1]
                pre = prefix.split('_PHOT')[0]
                cols = catalog.columns.str.rstrip('')
            all_cols = cols.to_list()
            fluxcols = all_cols[5:]
            flux_errs = []
            flux = []
            for f in fluxcols:
                if f.endswith('err'):
                    flux_errs.append(f)
                else:
                    flux.append(f)
            return path, pre, flux_errs, flux, new_ID
#data_array = 'CDFS247586'
#path, prefix, flux_errs, flux, new_ID = find_file('UDS_HST000149', 'fits')

#print(f'path:{ path}\nprefix:{ prefix}\n flux: {flux}\n flux_errs: {flux_errs}')

def load_vandels(object):
    path, prefix, flux_errs, flux_cols, new_ID = find_file(object, 'fits')
    cat_file = Table.read(path).to_pandas()
    catalog = pd.DataFrame(cat_file)
    ind = catalog.set_index(str(prefix)+ catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))# + catalog['CAT'].str.decode("utf-8"))

    fluxes = ind.loc[object, flux_cols].values
    fluxes = np.array(fluxes, dtype=float)
    fluxerrs = ind.loc[object,flux_errs].values
    fluxerrs = np.array(fluxerrs,dtype=float)
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


#print(load_vandels('CDFS_HST034930'))
