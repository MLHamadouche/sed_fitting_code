import numpy as np
import bagpipes as pipes
import pandas as pd
import LoadData as ld
from astropy.io import fits
from astropy.table import Table

exp = {}
exp["age"] = (0.1, 15.)
exp["tau"] = (0.3, 10.)
exp["massformed"] = (1., 15.)
exp["metallicity"] = (0., 2.5)

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0., 2.)

fit_instructions = {}
fit_instructions["redshift"] = (0., 10.)
fit_instructions["exponential"] = exp
fit_instructions["dust"] = dust

uds_hst_filt =  np.loadtxt("catalogs/UDS_HST_filt_list.txt", dtype="str")
uds_ground_filt = np.loadtxt("catalogs/UDS_GROUND_filt_list.txt", dtype = "str")
cdfs_ground_filt = np.loadtxt("catalogs/CDFS_GROUND_filt_list.txt", dtype="str")
cdfs_hst_filt= np.loadtxt("catalogs/CDFS_HST_filt_list.txt", dtype="str")

passive_cut = Table.read('FirstProjectCatalogs/xmatch_spec_derived237objs.fits').to_pandas()

objects = (passive_cut['FIELD'].str.decode("utf-8").str.rstrip() + passive_cut['ID_1'].astype(str).str.pad(6, side='left', fillchar='0')).to_list()
#print(objects)

filt_list = []
for object in objects:
    if 'CDFS-HST' in object:
        filt_list.append(cdfs_hst_filt)
    elif 'CDFS-GROUND' in object:
        filt_list.append(cdfs_ground_filt)
    elif 'UDS-HST' in object:
        filt_list.append(uds_hst_filt)
    else:
        filt_list.append(uds_ground_filt)

fit = pipes.fit_catalogue(objects,fit_instructions,ld.load_vandels, spectrum_exists=False,cat_filt_list=filt_list, make_plots=True, run="vandels_objects")
#
#fit_cat.plot_spectrum_posterior()  # Shows the input and fitted spectrum/photometry
#fit.plot_sfh_posterior()       # Shows the fitted star-formation history
#fit.plot_1d_posterior()        # Shows 1d posterior probability distributions
#fit.plot_corner()

fit.fit(verbose=False, mpi_serial=True)
