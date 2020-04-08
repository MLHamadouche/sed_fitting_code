import bagpipes as pipes
import numpy as np

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
