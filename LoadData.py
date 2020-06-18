import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import os
import re
from collections import OrderedDict
import adam_load_vandels as a_load

#fluxes and flux errors columns are in wavelength order same as filters
#must input ID as e.g. CDFS_GROUND0000xx where chars are filled to 6 padded by leading zeroes.

def find_file(ID, extension):
    new_ID = re.search('\d+', ID).group()
    new_ID = new_ID.lstrip('0')
    print(new_ID)
    for root, dirs, files in os.walk('/Users/PhDStuff/sed_fitting_code'):
        if "CDFS" in ID:
            #files = [filename for filename in files if filename.startswith("VANDELS_CDFS") and filename.endswith(".{ext}".format(ext=extension))]
            if 'HST' in ID:
                files = "VANDELS_CDFS_HST_PHOT_v1.0.fits"
            else:
                files = "VANDELS_CDFS_GROUND_PHOT_v1.0.fits"
        else:
            #files = [filename for filename in files if filename.startswith("VANDELS_UDS") and filename.endswith(".{ext}".format(ext=extension))]
            if 'HST' in ID:
                files = "VANDELS_UDS_HST_PHOT_v1.0.fits"
            else:
                files = "VANDELS_UDS_GROUND_PHOT_v1.0.fits"
        print(files)
        #for filename in files:
        path = os.path.join(root, files)
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
                    flux_errs = list(OrderedDict.fromkeys(flux_errs))

                    if f.endswith('_2as') or f.startswith('CH1') and f.endswith('tphot') or f.startswith('CH2') and f.endswith('tphot'):#f == 'isofactor'
                        flux.append(f)
                    flux = list(OrderedDict.fromkeys(flux))

            else:
                for fe in fluxcols:
                    if fe.endswith('_err'):
                        flux_errs.append(fe)
                    else:
                        flux.append(fe)

        return new_path, pre, flux_errs, flux, new_ID
#'CDFS_GROUND000013'
#path, prefix, flux_errs, flux, new_ID = find_file('CDFS-HST000013', 'fits')
#'CDFS_HST034930'
#print(f'path:{ path}\nprefix:{prefix}\n flux: {flux}\n flux_errs: {flux_errs}')

def load_vandels(object):
    path, prefix, flux_errs, flux_cols, new_ID = find_file(object, 'fits')
    print(path,prefix, flux_errs, flux_cols)
    cat_file = Table.read(path).to_pandas()
    catalog = pd.DataFrame(cat_file)
    #print(catalog['CAT'].str.decode('utf-8'))
    ind = catalog.set_index(str(prefix)+ catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))# + catalog['CAT'].str.decode("utf-8"))
    #print(ind)
    if 'isofactor' in catalog.columns:
        flux = []
        ferrs = []
        print('GROUND CATALOGUES')
        for f in flux_cols:
            iso = ind.loc[object, 'isofactor']
            if '_2as' in f:
                flux.append(ind.loc[object, f]*iso)
                offset = np.loadtxt("vandels/offsets_cdfs_ground.txt")
            else:
                flux.append(ind.loc[object, f])
                offset = np.loadtxt("vandels/offsets_uds_ground.txt")

        for fe in flux_errs:
            iso = ind.loc[object, 'isofactor']

            if '_2as' in fe:
                ferrs.append(ind.loc[object,fe]*iso)
                #offset = np.loadtxt("vandels/offsets_cdfs_ground.txt")
            else:
                ferrs.append(ind.loc[object,fe])


        photometry = np.c_[flux,ferrs]
        photometry[:,0]*=offset

    else:
        print('HST CATALOGUES')
        if 'UDS' in object:
            offset = np.loadtxt("vandels/offsets_uds_hst.txt")
        else:
            offset = np.loadtxt("vandels/offsets_cdfs_hst.txt")
        fluxes=ind.loc[object,flux_cols].values
        fluxerrs=ind.loc[object, flux_errs].values
        photometry = np.c_[fluxes,fluxerrs]
        #offsets = np.loadtxt("vandels/offsets_uds_hst.txt")
        photometry[:,0]*=offset

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
#print(load_vandels('UDS-HST035930'))
#print(a_load.load_vandels_phot('UDS-HST035930SELECT'))
#'CDFS_GROUND000013'
