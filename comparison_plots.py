import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import load_data as ld

print("{:e}".format(10**10.6))
input()
ross_objects = Table.read('/Users/massissiliahamadouche/Downloads/massi_cdfs_vandels_test_phot.fits').to_pandas()

objects = np.array('CDFS'+ ross_objects['ID'].astype(str).str.pad(6, side='left', fillchar='0'))
#for col in ross_objects.columns:
    #print(col)

my_results = Table.read('masstest2_ra_dec_catalogue.fits').to_pandas()

my_masses = my_results['mass']
ross_mass = ross_objects['log10(M*)']
massi_z = my_results['redshift']
ross_z =  ross_objects['ZSPEC']

matplotlib.rcParams['font.family'] = "AppleMyungjo"
#matplotlib.rcParams['font.family'] = "Hiragino Mincho Pro"
plt.scatter(ross_mass, np.log10(my_masses), color="red", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
plt.xlabel('Derived masses (10$^{m}$ solar masses)')
plt.ylabel("Massi's masses (10$^{m}$ solar masses)")
plt.savefig('MassComparisonPlot_testIGM.png')
plt.close()
#print(my_results)

plt.scatter(ross_z, massi_z, color="blue", s = 1)
plt.title(" Plot of Massi's results versus Ross's results")
plt.xlabel('Derived redshifts')
plt.ylabel("Massi's redshifts")
plt.savefig('RedshiftComparisonPlot_testIGM.png')
plt.close()
