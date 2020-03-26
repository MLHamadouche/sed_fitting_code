import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simps
from astropy.cosmology import FlatLambdaCDM
import spec_fluxes as sf

def photometry(waves, fluxes, redshift):

    filters = ["CH2", "HAWKI_K","ISAAC_Ks","CH1","VIMOS_U","f098m","f105w","f125w","f160w", "f435w","f606w", "f775w","f814w", "f850lp"]

    filter_curves = []

    for filter in filters:
        filter_curves.append(np.loadtxt("filters/"+str(filter)))
    #print(f'filter_curves.shape: {filter_curves.shape}')

    new_filter_curves = []
    #print(new_filter_curves.shape)
    eff_wavs = []
    print(f'redshift: {redshift.reshape(1,5)}')
    for f in filter_curves:
        flux_filter = f[:,1]/np.max(f[:,1])
        #print(f'flux_filter: {flux_filter.shape}')
        wav_filter = f[:,0]
        #print(f'wav_filter: {wav_filter.shape}')
        eff_wavs.append(np.sum(wav_filter*flux_filter)/np.sum(flux_filter))
        filter_interpol=(np.interp(waves*(1 + redshift.reshape(5,1)), wav_filter, flux_filter, left=0, right=0)) # shift Wavelengths by factor of redshift to the left to get lambda_obs
        #print(f'waves: {np.expand_dims(waves, axis =0)}')
        #print(f'waves * redshift: {np.interp(waves*(1 + redshift), wav_filter, flux_filter, left=0, right=0)}')
        #print(f'filter_interpol: {filter_interpol}')
        new_filter_curves.append(filter_interpol)

    print(f'new_filter_curves: {new_filter_curves}')
    c = 2.99792458*10**10 #cms/s cgs units
    #fluxes*=(3.826*10**33)#to get in ergs/s/Angstrom

    fluxes*=(3.826*10**33)
    print(fluxes.shape)
    # for conversion to ergs/s/cm^2/Angstrom

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    Mpc_m = 3.086*10**22 #Mpc in m
    Mpc_cm = Mpc_m*10**2 #in centimetres

    l_dist=cosmo.luminosity_distance(redshift).value*Mpc_cm #converting Mpc to cm
    lum_area = 4 *np.pi*((l_dist)**2)
    #print("l_dist=", l_dist)
    #print("lum area shape : ",lum_area.shape)
    f_lam = fluxes/(np.expand_dims(lum_area, axis=1))
    #print(f'flam.shape: {f_lam.shape}')
    #print(f_lam)

    new_fluxes = []
    # model photometry
    #print(f'new_fluxes.shape: {new_fluxes.shape}')
    waves=np.expand_dims(waves, axis = 0)#*np.expand_dims(np.ones(5), axis=1)
    #print(f'waves.shape:{waves.shape}')

    for i in range(len(redshift)):
        for m in range(14):
            #print(f'new_filter_curves[m].shape: {new_filter_curves[m].shape}')
            print(f'f_lam[i]:{f_lam[i]}')
            print(f'new_filter_curves[i,m]: {new_filter_curves[i][m].shape}')
            print(f'waves[i]: {waves[i]}')
            f_lam_x_filters = f_lam[i]*new_filter_curves[m]*waves[i]
            new_fluxes[m] =np.trapz(f_lam_x_filters, x=waves[i])/np.trapz(new_filter_curves[i][m]*waves[i], x=waves[i])
    #print(f'f_lam_x_filters: {f_lam_x_filters}')
    print(new_fluxes)

    return new_fluxes

#ages, waves,flux_grid = sf.spectrum()
#flux_grid = flux_grid[4][100]
#new_fluxes = photometry(waves, flux_grid, 1.04)
#print(new_fluxes*10**10)


#print(f'new_fluxes.shape: {new_fluxes}')

            #print(f'f_lam_x_filters: {f_lam_x_filters }')
            #print(f'f_lam_x_filters.shape: {f_lam_x_filters.shape}')
            #print(f'waves.shape:{waves[i].shape}')
            #print(f'new_filter_curves: {new_filter_curves.shape}')
            #print(f'(np.trapz(f_lam_x_filters, x=waves)):{(np.trapz(f_lam_x_filters, x=waves[i])).shape}')
            #print(f'(np.trapz(new_filter_curves[m]*waves, x=waves)): {(np.trapz(new_filter_curves[i][m]*waves[i], x=waves[i])).shape}')
