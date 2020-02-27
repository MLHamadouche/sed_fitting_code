import numpy as np
import matplotlib.pyplot as plt
import eff_wavs_filter_function as ewavs


# for the dust you need to just run over a grid of values for A(lambda)

#calzetti  law

def dust(eff_wavs):

    A_lam = np.zeros(eff_wavs.shape)

    for i in eff_wavs:
        if 6300 <= i <= 22000:
            A_lam = 4.05 + 2.659(-1.857 + (1.040/i))

        elif 1200 < i < 6300:
            A_lam = 4.05 + 2.659(-2.156 +1.509/i - 0.198/i**2 + 0.011/i**3)

    return A_lam
