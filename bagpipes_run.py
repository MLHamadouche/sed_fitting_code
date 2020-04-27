import bagpipes as pipes
import numpy as np
import load_data as ld
import eff_wavs_filter_function as ewavs
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import conversion_function as cf

burst = {}
burst["age"] = (0., 13.)           # Vary the age between 0 and 13 Gyr
burst["metallicity"] = 1.0  # Vary stellar metallicity between 0 and 5 times Solar
burst["massformed"] = (0., 13.)    # Vary the Log_10 of total mass formed from 0 to 13.

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0., 2.)

fit_instructions = {}
fit_instructions["redshift"] = (0., 10.)
fit_instructions["burst"] = burst
fit_instructions["dust"] = dust


filters = ["/Users/PhDStuff/sed_fitting_code/filters/CH2",
            "/Users/PhDStuff/sed_fitting_code/filters/HAWKI_K",
            "/Users/PhDStuff/sed_fitting_code/filters/ISAAC_Ks",
            "/Users/PhDStuff/sed_fitting_code/filters/CH1",
            "/Users/PhDStuff/sed_fitting_code/filters/VIMOS_U","/Users/PhDStuff/sed_fitting_code/filters/f098m",
            "/Users/PhDStuff/sed_fitting_code/filters/f105w","/Users/PhDStuff/sed_fitting_code/filters/f125w",
            "/Users/PhDStuff/sed_fitting_code/filters/f160w","/Users/PhDStuff/sed_fitting_code/filters/f435w","/Users/PhDStuff/sed_fitting_code/filters/f606w",
            "/Users/PhDStuff/sed_fitting_code/filters/f775w","/Users/PhDStuff/sed_fitting_code/filters/f814w", "/Users/PhDStuff/sed_fitting_code/filters/f850lp"]

#filters_list = np.loadtxt("//Users/massissiliahamadouche/anaconda3/lib/python3.7/site-packages/bagpipes/filters/massi.filt_list", dtype='str')

"""
model = pipes.model_galaxy(model_components, filt_list = filters)

#print(filters_list)

#model = pipes.model_galaxy(model_components, filt_list=filters_list)


fig = model.plot()
fig = model.sfh.plot()
"""

#model = pipes.model_galaxy(model_components, spec_wavs=np.arange(5000., 10000., 5.))

#fig = model.plot()

catalog = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

eff_wavs = ewavs.filter_wavs()

objects = np.array('CDFS'+ catalog['ID'].astype(str).str.pad(6, side='left', fillchar='0'))

flux_cols = ['CH2_flux', 'Ks_HAWKI_flux','Ks_ISAAC_flux','CH1_flux','VIMOS_U_flux','F098M_flux','F105W_flux','F125W_flux','F160W_flux','F435W_flux','F606W_flux','F775W_flux','F814W_flux','F850LP_flux']
flux_errs_cols = ['CH2_err', 'Ks_HAWKI_err','Ks_ISAAC_err','CH1_err','VIMOS_U_err','F098M_err','F105W_err','F125W_err','F160W_err', 'F435W_err','F606W_err', 'F775W_err','F814W_err', 'F850LP_err']


data = []

for i in range(len(objects)):
    data.append(objects[i])

data = data[0:5]

def load_data(data_array):
    # load up the relevant columns from the catalogue.
    # Find the correct row for the object we want.
    ind = catalog.set_index('CDFS' + catalog['ID'].astype(str).str.pad(6, side="left", fillchar="0"))
    # Extract the object we want from the catalogue.
    fluxes = ind.loc[data_array, flux_cols].values
    #print(cat_data)
    fluxerrs = ind.loc[data_array,flux_errs_cols].values
    print(fluxes)
    # Turn these into a 2D array.
    photometry = np.c_[fluxes, fluxerrs]

    # blow up the errors associated with any missing fluxes.
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
                              cat_filt_list=filters, run="guo_cat")

fit_cat.fit()
