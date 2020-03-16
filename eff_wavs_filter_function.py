import numpy as np
import matplotlib.pyplot as plt

def filter_wavs():

    filters = ["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]

    filter_curves = []
    for filter in filters:
        filter_curves.append(np.loadtxt("filters/"+str(filter)))

    new_filter_curves = []
        # need to shift spectrum by factor of redshift to get the overlap
    eff_wavs = []

    for f in filter_curves:
        flux_filter = f[:,1]/np.max(f[:,1])
        wav_filter = f[:,0]
        eff_wavs.append(np.sum(wav_filter*flux_filter)/np.sum(flux_filter))

    eff_wavs = np.array(eff_wavs)
    #print(eff_wavs)
    return eff_wavs
#eff_wavs = filter_wavs()
