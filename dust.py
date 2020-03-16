import numpy as np
import matplotlib.pyplot as plt
import eff_wavs_filter_function as ewavs
#calzetti dust model

<<<<<<< HEAD

=======
# for the dust you need to just run over a grid of values for A(lambda)

#calzetti  law
#e_wavs = ewavs.filter_wavs()
#wavs_mic = np.array(e_wavs)
#print(wavs_mic)
>>>>>>> 2abeb31ac48754858da9422362ec3673340cef52
def dust_masks(eff_wavs):

    mask0 = (eff_wavs < 1200)
    mask1 = (eff_wavs >= 6300) & (eff_wavs <= 31000)
    #wavs1 = wavs_mic[mask1]
<<<<<<< HEAD
    #print(f'wavs1:{wavs1}'
    mask2 = (eff_wavs < 6300) &  (eff_wavs >= 1200)

    wavs_mic = np.copy(eff_wavs)/10**4
=======
    #print(f'wavs1:{wavs1}')
    mask2 = (eff_wavs < 6300) &  (eff_wavs >= 1200)

    wavs_mic = np.copy(eff_wavs)/10**4

>>>>>>> 2abeb31ac48754858da9422362ec3673340cef52
    #wavs2 = wavs_mic[mask2]
    #print(f'wavs2:{wavs2}')

    k_lam = np.zeros_like(wavs_mic)
    #for i in wavs1:
    #    k_lam.append((4.05 + 2.659*(-1.857 + (1.040/i)))/4.05)
<<<<<<< HEAD

    k_lam[mask0] = ((wavs_mic[mask0]/0.12)**-0.77 * ((4.05 +2.659*(-2.156 + 1.509/0.12 - 0.198/0.12**2 + 0.011/0.12**3))))/4.05
    k_lam[mask1] = (4.05 + 2.659*(-1.857 + (1.040/wavs_mic[mask1])))/4.05
    k_lam[mask2] = (4.05 +2.659*(-2.156 + 1.509/wavs_mic[mask2] - 0.198/wavs_mic[mask2]**2 + 0.011/wavs_mic[mask2]**3))/4.05

    #plt.plot(wavs_mic,k_lam)
    #plt.xlabel("Wavelengths (\u03bc m)")
    #plt.ylabel(r'A$_\nu$/A$_\lambda$')
    #plt.xlim(0.1, 0.4)
    #plt.ylim(0., 5)
    #plt.show()
    return(k_lam)

#for j in wavs2:
#    k_lam.append((4.05 +2.659*(-2.156 + 1.059/i - 0.198/i**2 + 0.011/i**3))/4.05)

#print(f'k_lam2:{k_lam}')
=======
    k_lam[mask0] = ((wavs_mic[mask0]/0.12)**-0.77 * ((4.05 +2.659*(-2.156 + 1.509/0.12 - 0.198/0.12**2 + 0.011/0.12**3))))/4.05
    k_lam[mask1] = (4.05 + 2.659*(-1.857 + (1.040/wavs_mic[mask1])))/4.05
    k_lam[mask2] = (4.05 +2.659*(-2.156 + 1.509/wavs_mic[mask2] - 0.198/wavs_mic[mask2]**2 + 0.011/wavs_mic[mask2]**3))/4.05
    #print(f'k_lam:{k_lam}')

    #for j in wavs2:
    #    k_lam.append((4.05 +2.659*(-2.156 + 1.059/i - 0.198/i**2 + 0.011/i**3))/4.05)

    #print(f'k_lam2:{k_lam}')
    #plt.plot( wavs_mic,k_lam)
    #plt.xlim(0., 5.)
    #plt.show()
    return(k_lam)

#dust_masks(wavs_mic)
>>>>>>> 2abeb31ac48754858da9422362ec3673340cef52
