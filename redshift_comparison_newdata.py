import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from astropy.table import Table
from astropy.io import fits
from scipy import stats
import statistics


#all = Table.read('uds_validation_results_all.fits').to_pandas()
all = Table.read('ecdfs_validation_results_all.fits').to_pandas()
#ross_objects = Table.read('/Users/massi_zphot_test/ecdfs_zphot_validation_phot_zspec.fits').to_pandas()
#ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/UDS_zphot_validation_phot_fullinfo.fits').to_pandas()
#my_results = Table.read('uds_validation_catalogue.fits').to_pandas()
#my_results = Table.read('new_massi_UDS_validation_catalogue.fits').to_pandas()
my_results = Table.read('ECDFS_validation_catalogue.fits').to_pandas()
matplotlib.rcParams['font.family'] = "AppleMyungjo"

target = my_results['target']

massi_z = my_results['redshift']

spec_z =  all['zspec']
rjm_z = all['zphot_rjm']
barros = all['zphot_barros']
mortlock2 = all['zphot_alice_bc03']
ciras = all['zphot_ciras']
jarvis = all['zphot_jarvis']

def stats_z(z_phot, z_spec, N_obj):
    dz =(z_spec - z_phot)/(1+z_spec)
    counter = 0
    CO = abs(dz)
    no_CO = CO > 0.15
    Nco = np.sum(no_CO)

    CO_percentage = Nco/N_obj*100
    compute_MAD = stats.median_absolute_deviation(dz)
    sigma_dz_NMAD = compute_MAD
    print(np.sum(CO<=0.15))
    sigdznoCO=stats.tstd(dz[CO<=0.15])

    print(f'Nco: {Nco}, CO_percentage: {np.round(CO_percentage,3)}%,  sigma_dz: {np.round(sigma_dz_NMAD,5)}, sigdznoCO: {np.round(sigdznoCO,6)}')
    return Nco, CO_percentage,  dz, CO, sigma_dz_NMAD, sigdznoCO

Nco_m, CO_percentage_m,  dz_m, CO_m, sigma_dz_m,sigdz_m = stats_z(massi_z, spec_z, N_obj=1207)
Nco_c, CO_percentag_c,  dz_c, CO_c, sigma_dz_c, sigdz_c = stats_z(ciras, spec_z, N_obj=1207)
Nco, CO_percentage,  dz, CO, sigma_dz, sigdz = stats_z(mortlock2, spec_z, N_obj=1207)
Nco_r, CO_percentage_r, dz_r, CO_r, sigma_dz_r, sigdz_r = stats_z(rjm_z, spec_z, N_obj=1207)
Nco_b, CO_percentage_b, dz_b, CO_b, sigma_dz_b, sigdz_b = stats_z(barros, spec_z, N_obj=1207)
ncom, copercentm, dzm, com, sigdm, signom = stats_z(massi_z, mortlock2, N_obj = 1207)
ncoj, copercentj, dzj, coj, sigdj, signoj = stats_z(massi_z, jarvis, N_obj = 1207)
Nco_rm, CO_percentage_rm, dz_rm, CO_rm, sigma_dz_rm, sigdz_rm = stats_z(massi_z, rjm_z, N_obj=1207)

spc = np.random.uniform(0., 6.5, size=(1207))
#pht_plus = 0.85*spec_z - 0.15
#pht_neg  = 1.15*spec_z + 0.15
pht_plus = 0.85*rjm_z - 0.15
pht_neg  = 1.15*rjm_z+ 0.15
x=(pht_plus+0.15)/0.85
x2=(pht_neg-0.15)/1.15
spc1 = np.random.uniform(0., 6.5, size=(16))
spc2 = np.random.uniform(0., 6.5, size=(11))
#plt.scatter(ross_z, massi_z, color="m", s = 4)
plt.scatter(rjm_z, massi_z, color='darkseagreen', s = 12, marker='o', alpha = 0.9, edgecolors='black',  linewidth=0.5)
plt.scatter(rjm_z[dz_rm>0.15], massi_z[dz_rm>0.15], s=8, color ='lightcoral', marker='o',alpha=0.9, edgecolors= "black", linewidth=0.5)
plt.scatter(rjm_z[dz_rm<-0.15], massi_z[dz_rm<-0.15], s=8, color ='lightcoral', marker='o', alpha=0.9, edgecolors= "black", linewidth=0.5)
#plt.scatter(ross_z, outliers, color = 'k', s=4)
plt.plot(rjm_z, rjm_z, color="k", linewidth=0.7, label = '1:1 line')
plt.plot(rjm_z, pht_neg, ':', color = "k", linewidth=0.3)
plt.plot(rjm_z, pht_plus, ':', color = "k", linewidth=0.3)
#plt.xlim(0.,6.8)
#plt.ylim(0.,6.8)
plt.xlabel(r'rjm_z')
plt.ylabel(r'Massi_z')
props = dict(boxstyle='round', facecolor='darkseagreen', alpha=0.2)
#t = plt.text(0.1, 0.9, 'serif', size='small' )
textstr = '\n'.join((
    r'$N_{obj}$: $%s$' % (1207, ),
    r'$CO$%%: $%.3f$' % (CO_percentage_rm, ),
    r'$\sigma_{dz}: %.4f$' % sigma_dz_rm, ))

plt.text(0.2, 5., textstr, family='sans-serif', variant='normal', verticalalignment='top', bbox=props)
plt.title(r'ECDFS Validation Comparison')
plt.savefig('ECDFS_rjm_validation_redshifts.png')
plt.close()
