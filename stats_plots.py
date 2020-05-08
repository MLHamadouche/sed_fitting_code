import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import bagpipes as pipes
from scipy import stats

ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()
bagpipes = Table.read('pipes/cats/guo_cat.fits').to_pandas()
objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
#my_results = Table.read('full_igm_ra_dec_catalogue.fits').to_pandas()
my_results = Table.read('stellar_mass_ra_dec_catalogue.fits').to_pandas()

matplotlib.rcParams['font.family'] = "AppleMyungjo"
########### AGES ####################################

massi_ages = my_results['age']/10**9
bagpipes_ages = bagpipes['burst:age_50']

###### targets  #####################################
ID_pipes  = bagpipes['#ID']
target = my_results['target']
########### DUST ####################################
my_dust = my_results['dust']
ross_dust = ross_objects['AV']
bagpipes_dust = bagpipes['dust:Av_50']
####### MASSES ######################################
#my_masses = np.log10(my_results['mass'])
massi_mass = np.log10(my_results['stellar_mass'])
massi_formedmass = np.log10(my_results['formed_mass'])
ross_mass = ross_objects['log10(M*)']
bagpipes_mass = bagpipes['burst:massformed_50']
delta_mass = ross_mass - massi_mass
delta_bp_mass = bagpipes_mass - massi_mass
bagpipes_stellar_mass = bagpipes['stellar_mass_50']

####### REDSHIFT ####################################
massi_z = my_results['redshift']
bagpipes_z = bagpipes['redshift_50']
ross_z =  ross_objects['ZSPEC']

def stats_z(my_z, spec_z, N_obj):
    z_spec1 = spec_z
    z_phot = my_z
    dz = (z_spec1 - z_phot)/(1+z_spec1)
    counter = 0
    CO = abs(dz)
    no_CO = []
    for i in range(N_obj):
        if CO[i] > 0.15:
            counter+=1
            no_CO.append(counter)

        Nco = len(no_CO)

        CO_percentage = Nco/N_obj *100

    compute_MAD = stats.median_absolute_deviation(dz)
    print(f'Nco: {Nco}, CO_percentage: {np.round(CO_percentage,3)}%, compute_MAD: {np.round(compute_MAD,5)}')
    return Nco, CO_percentage, compute_MAD, dz

Nco, CO_percentage, compute_MAD, dz = stats_z(massi_z, ross_z, N_obj=309)
Ncob, CO_percentageb, compute_MADb, dzb = stats_z(massi_z, bagpipes_z, N_obj=309)

spc = np.random.uniform(0.1, 7., size=(309))
pht_plus = 0.85*ross_z - 0.15#randint(0., 7., 309)
pht_neg  = 1.15*ross_z - 0.15
pht_plus2 = 0.85*bagpipes_z - 0.15#randint(0., 7., 309)
pht_neg2  = 1.15*bagpipes_z - 0.15
pht_plus3 = 0.85*ross_z - 0.15#randint(0., 7., 309)
pht_neg3  = 1.15*ross_z - 0.15
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
ax1.scatter(ross_z, massi_z, color="m", s = 2)
ax2.scatter(bagpipes_z, massi_z, color="m", s = 2)
ax3.scatter(ross_z, bagpipes_z, color="m", s = 2)
l1 = ax1.plot(ross_z, ross_z, color="k", linewidth=0.7, label = '1:1 line')
l2 = ax2.plot(bagpipes_z, bagpipes_z,color="k", linewidth=0.7, label = '1:1 line' )
l3 = ax3.plot(ross_z, ross_z,color="k", linewidth=0.7, label = '1:1 line' )
l4 = ax1.plot(ross_z, pht_plus, '--', color = "k", linewidth=0.7)
l5 = ax1.plot(ross_z, pht_neg, '--', color = "k", linewidth=0.7)
l4 = ax2.plot(bagpipes_z, pht_plus2, '--', color = "k", linewidth=0.7)
l5 = ax2.plot(bagpipes_z, pht_neg2, '--', color = "k", linewidth=0.7)
l4 = ax3.plot(ross_z, pht_plus3, '--', color = "k", linewidth=0.7)
l5 = ax3.plot(ross_z, pht_neg3, '--', color = "k", linewidth=0.7)
#plt.title(" Plot of Massi's results versus cat results")
ax1.set_xlabel('Ross z')
ax2.set_ylabel("Massi's redshifts")
ax3.set_xlabel('Ross z')
ax3.set_ylabel("Bagpipes z")
ax2.set_xlabel('Bagpipes z')
ax1.set_ylabel("Massi's redshifts")
fig.suptitle("Plot of redshift results")
plt.savefig('_all_z_results.png')
plt.close()

############## dust v dust #########################################
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,6))
ax1.scatter(ross_dust, my_dust, color="m", s = 2)
ax2.scatter(bagpipes_dust, my_dust, color="m", s = 2)
l1 = ax1.plot(ross_dust, ross_dust, color="k", linewidth=0.7, label = '1:1 line')
l2 = ax2.plot(bagpipes_dust, bagpipes_dust,color="k", linewidth=0.7, label = '1:1 line' )
#plt.title(" Plot of Massi's results versus cat results")
ax1.set_xlabel('Ross dust (A_v)')
ax2.set_ylabel("Massi's dust (A_v) ")
ax2.set_xlabel('Bagpipes dust (A_v)')
ax1.set_ylabel("Massi's dust (A_v)")
fig.suptitle("Plot of Massi's results versus bagpipes & Ross' results")
plt.savefig('my_res_v_b&r_dust_results.png')
plt.close()

############## age v age with mass cbar ###############################
fig, ax = plt.subplots()
c = massi_mass # color of points
im = ax.scatter(bagpipes_ages, massi_ages,c=c, s=5, linewidth=0.7, cmap=plt.cm.BuPu_r)
im2 = ax.plot(bagpipes_ages, bagpipes_ages,linewidth=0.7, color="black")
# Add a colorbar
plt.xlabel(r'bagpipes_age 10$^{9}$ yrs')
plt.ylabel(r'massi_age 10$^{9}$ yrs')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'Mass(10$^m$$M_{\odot}$)')
plt.xlim(0,1.5)
plt.ylim(0,2)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('age_mass_log10_cbar.png')
plt.close()

######### mass v mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(bagpipes_mass, massi_mass,c=c,s=5, cmap=plt.cm.BuPu_r)
im2 = ax.plot(bagpipes_stellar_mass, bagpipes_stellar_mass, linewidth=0.7, color="black")
plt.xlabel(r'bagpipes_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.xlim(8,11.5)
plt.ylim(7.5,11.5)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('pipesmassvmass_age_cbar.png')
plt.close()

######### mass v mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(ross_mass, massi_mass, c=c,s=5, cmap=plt.cm.BuPu_r)
im2 = ax.plot(ross_mass, ross_mass,linewidth=0.7, color="black")
plt.xlabel(r'ross_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.title(" Plot of Massi's results versus Ross results")
plt.savefig('stellarrossmassvmass_age_cbar.png')
plt.close()

######### mass v stellar-mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(bagpipes_stellar_mass, massi_mass,c=c,s=5, cmap=plt.cm.BuPu_r)
im2 = ax.plot(bagpipes_stellar_mass, bagpipes_stellar_mass, linewidth=0.7, color="black")
plt.xlabel(r'bagpipes_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.xlim(8,11.5)
plt.ylim(7.5,11.5)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('pipesstellarmassvmass_age_cbar.png')
plt.close()
