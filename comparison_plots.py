import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import load_data as ld
import scipy.optimize as opt

print("{:e}".format(10**10.6))
input()
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
#for col in ross_objects.columns:
    #print(col)

my_results = Table.read('full_igm_ra_dec_catalogue.fits').to_pandas()

my_masses = my_results['mass']
ross_mass = ross_objects['log10(M*)']
massi_z = my_results['redshift']
ross_z =  ross_objects['ZSPEC']
delta_mass = ross_mass - np.log10(my_masses)

matplotlib.rcParams['font.family'] = "AppleMyungjo"
#matplotlib.rcParams['font.family'] = "Hiragino Mincho Pro"
plt.scatter(ross_mass, np.log10(my_masses), color="red", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
plt.xlabel('Derived masses (10$^{m}$ solar masses)')
plt.ylabel("Massi's masses (10$^{m}$ solar masses)")
plt.savefig('MassComparisonPlot_IGM.png')
plt.close()
#print(my_results)

plt.scatter(ross_mass,delta_mass,  color="red", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
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

best_chisq = opt.minimize(chisq,x0,args=[ross_z, massi_z, yerrs])

print(best_chisq)
m = 1.01
c = 1.4901161193847656e-08
plt.plot(ross_z, model(m,c, ross_z),'',label='model')
#plt.errorbar(ross_z, massi_z,yerr=yerrs,label='errors', ls=" ")
plt.scatter(ross_z, massi_z, color="blue", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
plt.xlabel('Derived redshifts')
plt.ylabel("Massi's redshifts")
plt.savefig('RedshiftComparisonPlottest_IGM.png')
plt.close()
