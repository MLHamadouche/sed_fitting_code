import bagpipes as pipes
import numpy as np
import load_data as ld
import eff_wavs_filter_function as ewavs


exp = {}                          # Tau model star formation history component
exp["age"] = 2.                   # Gyr
exp["tau"] = 0.75                 # Gyr
exp["massformed"] = 9.            # log_10(M*/M_solar)
exp["metallicity"] = 0.5          # Z/Z_oldsolar

dust = {}                         # Dust component
dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
dust["Av"] = 0.4                  # magnitudes

model_components = {}                   # The model components dictionary
model_components["redshift"] = 1.13     # Observed redshift
model_components["exponential"] = exp
model_components["dust"] = dust


filters = ["filters/filters_m/CH2",
            "filters/filters_m/HAWKI_K",
            "filters/filters_m/ISAAC_Ks",
            "filters/filters_m/CH1",
            "filters/filters_m/VIMOS_U","filters/filters_m/f098m",
            "filters/filters_m/f105w","filters/filters_m/f125w",
            "filters/filters_m/f160w","filters/filters_m/f435w","filters/filters_m/f606w",
            "filters/filters_m/f775w","filters/filters_m/f814w", "filters/filters_m/f850lp"]

filters_list = np.loadtxt("//Users/massissiliahamadouche/anaconda3/lib/python3.7/site-packages/bagpipes/filters/massi.filt_list", dtype='str')


model = pipes.model_galaxy(model_components, filt_list = filters)

#print(filters_list)

model = pipes.model_galaxy(model_components, filt_list=filters_list)


fig = model.plot()
fig = model.sfh.plot()


#model = pipes.model_galaxy(model_components, spec_wavs=np.arange(5000., 10000., 5.))

#fig = model.plot()

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


eff_wavs = ewavs.filter_wavs()
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))

data = []
for i in range(len(objects)):
    data.append(objects[i])
#print(data)
ID, fluxes_obs_raw, fluxerrs_obs_raw = ld.load_catalog_data(data)

fluxes_obs = cf.conversion_func(fluxes_obs_raw, eff_wavs)
fluxerrs_obs = cf.conversion_func(fluxerrs_obs_raw, eff_wavs)


fit_cat = pipes.fit_catalogue(IDs, fit_instructions, load_goodss, spectrum_exists=False,
                              cat_filt_list=filters, run="guo_cat")

fit_cat.fit(verbose=False)
