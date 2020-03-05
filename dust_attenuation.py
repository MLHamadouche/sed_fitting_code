import numpy as np
import matplotlib.pyplot as plt
import eff_wavs_filter_function as ewavs


# for the dust you need to just run over a grid of values for A(lambda)

#calzetti  law

#print(eff_wavs)
def dust(eff_wavs):
    k_lam =[]

    for i in eff_wavs:

        k_lam.append((4.05 + 2.659*(-1.857 + (1.040/i)))/4.05)

    return(k_lam)
#e_wavs = ewavs.filter_wavs()
#eff_wavs = np.array(e_wavs)
#k_lam_ = dust(eff_wavs)
#print(k_lam_)
