import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import os
import re
from collections import OrderedDict

#fluxes and flux errors columns are in wavelength order same as filters
#must input ID as e.g. CDFS_GROUND0000xx where chars are filled to 6 padded by leading zeroes.

def find_file(ID, extension):
    new_ID = re.search('\d+', ID).group()
    new_ID = new_ID.lstrip('0')
    print(new_ID)
    for root, dirs, files in os.walk('/Users/PhDStuff/sed_fitting_code'):
        if ID.startswith("CDFS"):
            files = [filename for filename in files if filename.startswith("VANDELS_CDFS") and filename.endswith(".{ext}".format(ext=extension))]
        else:
            files = [filename for filename in files if filename.startswith("VANDELS_UDS") and filename.endswith(".{ext}".format(ext=extension))]

        for filename in files:
            path = os.path.join(root, filename)
            #print(path)
            data_file = Table.read(path).to_pandas()
            catalog = pd.DataFrame(data_file)

            ID_list = catalog['ID'].astype(str)

            if new_ID in ID_list.values:
                new_path = path
                both = os.path.split(path)
                prefix_for_load = both[1]
                prefix = prefix_for_load.split('VANDELS_')[1]
                pre = prefix.split('_PHOT')[0]
                cols = catalog.columns.str.rstrip('')
                all_cols = cols.to_list()
                fluxcols = all_cols[5:]
                flux_errs = []
                flux = []
                if 'isofactor' in catalog.columns:
                    for f in fluxcols:
                        if f.endswith('_2as_err') or f.startswith('CH1') and f.endswith('tphot_err') or f.startswith('CH2') and f.endswith('tphot_err'): #f == 'isofactor' and
                            flux_errs.append(f)
                        #else:
                            #flux_errs.remove(f)
                        flux_errs = list(OrderedDict.fromkeys(flux_errs))

                        if f.endswith('_2as') or f.startswith('CH1') and f.endswith('tphot') or f.startswith('CH2') and f.endswith('tphot'):#f == 'isofactor'
                            flux.append(f)
                        #else:
                            #flux.remove(f)
                        flux = list(OrderedDict.fromkeys(flux))

                else:
                    for fe in fluxcols:
                        if fe.endswith('_err'):
                            flux_errs.append(fe)
                        else:
                            flux.append(fe)

        return new_path, pre, flux_errs, flux, new_ID
#'CDFS_GROUND000013'
#path, prefix, flux_errs, flux, new_ID = find_file('CDFS_GROUND000013', 'fits')
#'CDFS_HST034930'
#print(f'path:{ path}\nprefix:{prefix}\n flux: {flux}\n flux_errs: {flux_errs}')
def load_vandels(object):
    path, prefix, flux_errs, flux_cols, new_ID = find_file(object, 'fits')
    print(path,prefix, flux_errs, flux_cols)
    cat_file = Table.read(path).to_pandas()
    catalog = pd.DataFrame(cat_file)
    ind = catalog.set_index(str(prefix)+ catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))# + catalog['CAT'].str.decode("utf-8"))

    if 'isofactor' in catalog.columns:
        flux = []
        ferrs = []
        print('GROUND CATALOGUES')
        for f in flux_cols:
            iso = ind.loc[object, 'isofactor']

            if '_2as' in f:
                flux.append(ind.loc[object, f]*iso)
            else:
                flux.append(ind.loc[object, f])

        for fe in flux_errs:
            iso = ind.loc[object, 'isofactor']

            if '_2as' in fe:
                ferrs.append(ind.loc[object,fe]*iso)
            else:
                ferrs.append(ind.loc[object,fe])

        photometry = np.c_[flux,ferrs]

    else:
        print('HST CATALOGUES')

        fluxes=ind.loc[object,flux_cols].values
        fluxerrs=ind.loc[object, flux_errs].values
        photometry = np.c_[fluxes,fluxerrs]

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

#ID = 'UDS_HST035930'
#print(load_vandels(ID))
#'CDFS_GROUND000013'
