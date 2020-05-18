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
#my_results = Table.read('stellar_mass_ra_dec_catalogue.fits').to_pandas()
my_results = Table.read('new_mass_fixed_catalogue.fits').to_pandas()
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
massi_stel_mass = np.log10(my_results['stellar mass'])
massi_formedmass = np.log10(my_results['formed mass'])
ross_mass = ross_objects['log10(M*)']
bagpipes_formed_mass = bagpipes['burst:massformed_50']
delta_mass = ross_mass - massi_stel_mass
bagpipes_stellar_mass = bagpipes['stellar_mass_50']
delta_bp_mass = bagpipes_stellar_mass - massi_stel_mass
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
    sigma_dz = compute_MAD/1.4286
    print(f'Nco: {Nco}, CO_percentage: {np.round(CO_percentage,3)}%, compute_MAD: {np.round(compute_MAD,5)},  compute_MAD: {np.round(sigma_dz,5)}')
    return Nco, CO_percentage, compute_MAD, dz, sigma_dz

Nco, CO_percentage, compute_MAD, dz, sig_dz = stats_z(massi_z, ross_z, N_obj=309)
Ncob, CO_percentageb, compute_MADb, dzb, sig_dzb = stats_z(massi_z, bagpipes_z, N_obj=309)
Ncobr, CO_percentagebr, compute_MADbr, dzbr, sig_dzrb = stats_z(bagpipes_z, ross_z, N_obj=309)

spc = np.random.uniform(0., 6., size=(309))
pht_plus = 0.85*spc - 0.15#randint(0., 7., 309)
pht_neg  = 1.15*spc +0.15
pht_plus2 = 0.85*spc- 0.15#randint(0., 7., 309)
pht_neg2  = 1.15*spc + 0.15
pht_plus3 = 0.85*spc - 0.15#randint(0., 7., 309)
pht_neg3  = 1.15*spc + 0.15
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
ax1.scatter(ross_z, massi_z, color='mediumorchid', s = 10, marker='o', alpha=0.3)
ax1.scatter(ross_z[dz>0.15], massi_z[dz>0.15], s=10, color ='rebeccapurple', marker='o',)
ax1.scatter(ross_z[dz<-0.15], massi_z[dz<-0.15], s=10, color ='rebeccapurple', marker='o')
ax2.scatter(bagpipes_z, massi_z, color='mediumorchid', s = 10, marker='o', alpha=0.3)
ax2.scatter(bagpipes_z[dzb>0.15], massi_z[dzb>0.15], s=10, color ='rebeccapurple', marker='o',)
ax2.scatter(bagpipes_z[dzb<-0.15], massi_z[dzb<-0.15], s=10, color ='rebeccapurple', marker='o')
ax3.scatter(ross_z, bagpipes_z, color='mediumorchid', s = 10, marker='o', alpha=0.3)
ax3.scatter(ross_z[dzbr>0.15], bagpipes_z[dzbr>0.15], s=10, color ='rebeccapurple', marker='o',)
ax3.scatter(ross_z[dzbr<-0.15], bagpipes_z[dzbr<-0.15], s=10, color ='rebeccapurple', marker='o')
l1 = ax1.plot(spc, spc, color="k", linewidth=0.8, label = '1:1 line')
l2 = ax2.plot(spc, spc,color="k", linewidth=0.8, label = '1:1 line' )
l3 = ax3.plot(spc, spc,color="k", linewidth=0.8, label = '1:1 line' )
l4 = ax1.plot(spc, pht_plus, '--', color = "k", linewidth=0.5)
l5 = ax1.plot(spc, pht_neg, '--', color = "k", linewidth=0.5)
l4 = ax2.plot(spc, pht_plus2, '--', color = "k", linewidth=0.5)
l5 = ax2.plot(spc, pht_neg2, '--', color = "k", linewidth=0.5)
l4 = ax3.plot(spc, pht_plus3, '--', color = "k", linewidth=0.5)
l5 = ax3.plot(spc, pht_neg3, '--', color = "k", linewidth=0.5)
#plt.title(" Plot of Massi's results versus cat results")
ax1.set_xlabel('Ross z')
ax2.set_ylabel("Massi's redshifts")
ax3.set_xlabel('Ross z')
ax3.set_ylabel("Bagpipes z")
ax2.set_xlabel('Bagpipes z')
ax1.set_ylabel("Massi's redshifts")
props = dict(boxstyle='round', facecolor='plum', alpha=0.4)
textstr1 = '\n'.join((
    r'N$_{obj}$: $%s$' % (309, ),
    r'$\mathrm{CO}: %.3f$%%' % (CO_percentage, ),
    r'$\sigma_{dz}: %.4f$' % (sig_dz, )))
ax1.text(0.2, 7.,textstr1, family='serif',variant ='normal',size='x-small', verticalalignment='top', bbox=props)
textstr2 = '\n'.join((
    r'N$_{obj}$: $%s$' % (309, ),
    r'$\mathrm{CO}: %.3f$%%' % (CO_percentageb, ),
    r'$\sigma_{dz}: %.4f$' % (sig_dzb, )))
ax2.text(0.2, 7.,textstr2, family='serif',variant ='normal',size='x-small', verticalalignment='top', bbox=props)
textstr3 = '\n'.join((
    r'N$_{obj}$: $%s$' % (309, ),
    r'$\mathrm{CO}: %.3f$%%' % (CO_percentagebr, ),
    r'$\sigma_{dz}: %.4f$' % (sig_dzrb, )))
ax3.text(0.2, 7.,textstr3, family='serif',variant ='normal', size='x-small', verticalalignment='top', bbox=props)
fig.suptitle("Plot of redshift results")
plt.savefig('new_all_z_results.png')
plt.close()
#plt.show()
############## dust v dust #########################################
dust1to1 = np.random.uniform(0., 2., size=(309))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,6))
ax1.scatter(ross_dust, my_dust, color='mediumorchid', s = 15, alpha=0.7)
ax2.scatter(bagpipes_dust, my_dust, color='mediumorchid', s = 15, alpha=0.7)
l1 = ax1.plot(dust1to1, dust1to1, color="k", linewidth=0.7, label = '1:1 line')
l2 = ax2.plot(dust1to1, dust1to1,color="k", linewidth=0.7, label = '1:1 line' )
#plt.title(" Plot of Massi's results versus cat results")
ax1.set_xlabel('Ross dust (A_v)')
ax2.set_ylabel("Massi's dust (A_v) ")
ax2.set_xlabel('Bagpipes dust (A_v)')
ax1.set_ylabel("Massi's dust (A_v)")
fig.suptitle("Plot of Massi's results versus bagpipes & Ross' results")
plt.savefig('new_dust_results.png')
plt.close()
#plt.show()
############## age v age with mass cbar ###############################
massi_age = my_results['age']
bagpipes_age = bagpipes['burst:age_50']*10**9
fig, ax = plt.subplots()
c = massi_stel_mass # color of points
im = ax.scatter(bagpipes_age, massi_age,c=c, s=15, linewidth=0.3, cmap=plt.cm.BuPu_r,  alpha=0.7)
im2 = ax.plot(bagpipes_age, bagpipes_age,linewidth=0.3, color="black")
# Add a colorbar
plt.xlabel(r'bagpipes_age 10$^{9}$ yrs')
plt.ylabel(r'massi_age 10$^{9}$ yrs')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'Mass(10$^m$$M_{\odot}$)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10**7, 10**9.5)
plt.ylim(10**7.5, 10**9.5)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('new_age_mass_log10_cbar.png')
plt.close()
#plt.show()
######### mass v mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(bagpipes_formed_mass, massi_formedmass,c=c,s=15, cmap=plt.cm.BuPu_r,  alpha=0.7)
im2 = ax.plot(bagpipes_formed_mass, bagpipes_formed_mass, linewidth=0.7, color="black")
plt.xlabel(r'bagpipes_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_formed_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.xlim(8.,11.5)
plt.ylim(8.,11.5)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('new_pipes_formedmassv.png')
plt.close()
#plt.show()
######### mass v mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(ross_mass, massi_stel_mass, c=c,s=15, cmap=plt.cm.BuPu_r,  alpha=0.7)
im2 = ax.plot(ross_mass, ross_mass,linewidth=0.7, color="black")
plt.xlabel(r'ross_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_stel_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.title(" Plot of Massi's results versus Ross results")
plt.savefig('new_rossmassvmass_age_cbar.png')
plt.close()
#plt.show()
######### mass v stellar-mass with age cbar ####################################
fig, ax = plt.subplots()
c = massi_ages*10**9  # color of points
im = ax.scatter(bagpipes_stellar_mass, massi_stel_mass,c=c,s=15, cmap=plt.cm.BuPu_r,  alpha=0.7)
im2 = ax.plot(bagpipes_stellar_mass, bagpipes_stellar_mass, linewidth=0.7, color="black")
plt.xlabel(r'bagpipes_stel_mass 10$^{m}$ $M_{\odot}$')
plt.ylabel(r'massi_stel_mass 10$^{m}$ $M_{\odot}$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'age(10$^9$yrs)')
plt.xlim(8.,11.5)
plt.ylim(8.,11.5)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.savefig('new_pipes_stellarmass_cbar.png')
plt.close()
#plt.show()
