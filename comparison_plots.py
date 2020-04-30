import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import load_data as ld
import scipy.optimize as opt

#print("{:e}".format(10**10.6))

ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()
bagpipes = Table.read('pipes/cats/guo_cat.fits').to_pandas()
objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
#for col in ross_objects.columns:
    #print(col)
pipes  = bagpipes['#ID']

my_results = Table.read('full_igm_ra_dec_catalogue.fits').to_pandas()
my_dust = my_results['dust']
bagpipes_dust = bagpipes['dust:Av_50']
my_masses = np.log10(my_results['mass'])
ross_mass = ross_objects['log10(M*)']
massi_z = my_results['redshift']
bagpipes_mass = bagpipes['burst:massformed_50']
bagpipes_z = bagpipes['redshift_50']
ross_z =  ross_objects['ZSPEC']
delta_mass = ross_mass - my_masses
delta_bp_mass = bagpipes_mass - my_masses
matplotlib.rcParams['font.family'] = "AppleMyungjo"
#matplotlib.rcParams['font.family'] = "Hiragino Mincho Pro"
plt.scatter(bagpipes_dust, my_dust,  color="blue", s = 1)
plt.title(" Plot of Massi's results versus Bagpipes results")
plt.xlabel('Bagpipes dust (A_v)')
plt.ylabel("Massi's dust (A_v)")
plt.savefig('dust_bp_comparison.png')
plt.close()
#print(my_results)
"""
plt.scatter(ross_mass,delta_mass,  color="red", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
plt.xlabel('Derived masses (10$^{m}$ solar masses)')
plt.ylabel("Delta mass (10$^{m}$ solar masses)")
plt.savefig('DeltaMassComparisonPlot_IGM.png')
plt.close()
"""

plt.scatter(bagpipes_mass, delta_mass,  color="red", s = 1)
plt.title(" Plot of Massi's results versus Bagpipes results")
plt.xlabel('Derived masses (10$^{m}$ solar masses)')
plt.ylabel("Delta mass (10$^{m}$ solar masses)")
plt.savefig('DeltaMassComparisonPlot_IGM.png')
plt.close()

def model(a,b, x1):
    return a*x1 + b

def chisq(param, args):

    m, c = param
    x1, y1, yerrs = args
    y_model = model(m,c, x1)

    chi2 = np.sum(((y1 - y_model)/yerrs)**2)
    print(m, c, chi2)
    return chi2

yerrs = np.ones(309)

x0 = [0., 0.]

print(x0)
#CDFS034186
#CDFS019527
print(np.max(bagpipes_mass))

best_chisq_m = opt.minimize(chisq,[6.5,7],args=[bagpipes_mass, my_masses, yerrs])

mm, mc = best_chisq_m.x
print(f'mm:{mm}, mc: {mc}')
plt.plot(bagpipes_mass, model(mm, mc, bagpipes_mass),'',label='model')
#(mm:0.7924392193530521, mc: 1.5071569509261304), chisq: 35.997870567558394
#plt.plot(bagpipes_mass, model(mm, mc, bagpipes_mass),'',label='model')
#plt.scatter(bagpipes_mass, my_masses, color="red", s = 1)
plt.scatter(bagpipes_mass, my_masses, color="red", s = 1)
#plt.plot(bagpipes_mass, model(mm, mc, bagpipes_mass),'',label='model')
#plt.scatter(bagpipes_mass, my_masses, color="red", s = 1)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.xlabel('Bagpipes masses (10$^{m}$ solar masses)')
plt.ylabel("Massi's masses (10$^{m}$ solar masses)")

#plt.show()
plt.savefig('MassComparisonPlotbagpipes_IGM.png')
plt.close()


for i in range(len(bagpipes_z)):
    if bagpipes_z[i] > 1:
        best_chisq_z = opt.minimize(chisq,x0,args=[bagpipes_z[i], massi_z, yerrs])

zm, zc = best_chisq_z.x

#plt.errorbar(ross_z, massi_z,yerr=yerrs,label='errors', ls=" ")

plt.plot(bagpipes_z, model(zm, zc, bagpipes_z),'', color="red", label='model')
plt.scatter(bagpipes_z, massi_z, color="blue", s = 1)
plt.title(" Plot of Massi's results versus bagpipes results")
plt.xlabel('Bagpipes redshifts')
plt.ylabel("Massi's redshifts")
plt.savefig('RedshiftComparison_bagpipes_IGM.png')
plt.close()

best_chisq_rossz = opt.minimize(chisq,x0,args=[ross_z, massi_z, yerrs])

mr, cr = best_chisq_rossz.x
plt.plot(ross_z, model(mr, cr, ross_z),'', color="red", label='model')
plt.scatter(ross_z, massi_z, color="blue", s = 1)
plt.title(" Plot of Massi's results versus Ross' results")
plt.xlabel('ZSPEC')
plt.ylabel("Massi's redshifts")
plt.savefig('RedshiftComparison_Ross_IGM.png')
plt.close()

best_chisq_rossm = opt.minimize(chisq,[1.,1.],args=[ross_mass, my_masses, yerrs])

rm, rc= best_chisq_rossm.x
m, c = 1.01, 1.4901161193847656*10**-8
chis = 37.66017750019812
plt.plot(ross_mass, model(rm, rc, ross_mass),'', color="red", label='model')
plt.scatter(ross_mass, my_masses, color="blue", s = 1)
plt.title(" Plot of Massi's results versus Ross' results")
plt.xlabel("Ross' masses")
plt.ylabel("Massi's masses")
plt.savefig('MassComparison_Ross_IGM.png')
plt.close()
