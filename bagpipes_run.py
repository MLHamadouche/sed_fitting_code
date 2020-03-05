import bagpipes as pipes
import numpy as np

exp = {}                          # Tau model star formation history component
exp["age"] = 3.                   # Gyr
exp["tau"] = 0.75                 # Gyr
exp["massformed"] = 9.            # log_10(M*/M_solar)
exp["metallicity"] = 0.5          # Z/Z_oldsolar

dust = {}                         # Dust component
dust["type"] = "Calzetti"         # Define the shape of the attenuation curve
dust["Av"] = 0.2                  # magnitudes

model_components = {}                   # The model components dictionary
model_components["redshift"] = 1.0      # Observed redshift
model_components["exponential"] = exp
model_components["dust"] = dust


#/Users/massissiliahamadouche/anaconda3/lib/python3.7/site-packages/bagpipes/filters
#filters = ["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]

filters_list = np.loadtxt("//Users/massissiliahamadouche/anaconda3/lib/python3.7/site-packages/bagpipes/filters/massi.filt_list", dtype='str')

#print(filters_list)

model = pipes.model_galaxy(model_components, filt_list=filters_list)

fig = model.plot()
fig = model.sfh.plot()


model = pipes.model_galaxy(model_components, filt_list=filters_list, spec_wavs=np.arange(5000., 10000., 5.))

fig = model.plot()
