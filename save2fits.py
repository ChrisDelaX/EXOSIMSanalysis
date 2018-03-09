import numpy as np
from astropy.io import fits

ar = np.loadtxt('/Users/cdelacroix/Desktop/Bijan/QEtable.txt')

names = ['Basic_Multi2_DD', 'Astro_Multi2_DD', 'Astro_Multi2_std_Si', 
        'Basic_Multi2_std_Si', 'Astro_Midband_DD', 'Basic_NIR_DD', 'Basic_Midband_Std_Si']

for i,name in enumerate(names):
    dat = np.array([ar[:,0],ar[:,i+1]])
    hdu = fits.PrimaryHDU(dat)
    fits.HDUList([hdu]).writeto(name+'.fits')

###

import numpy as np
from astropy.io import fits

perfs = ['mean_intensity', 'contrast', 'thruput', 'area', 'occ_trans']
cols = [2, 3, 4, 6, 7]
path = '/Users/cdelacroix/Desktop/Bijan/'

#name = 'B25_FIT_770'#'B25_FIT_660'#'L3_FIT_770'#'L3_FIT_660'
name = 'B22_FIT_565'#'G22_FIT_565'

ar = np.loadtxt(path+name+'.txt')

for i, perf in enumerate(perfs):
    hdu = fits.PrimaryHDU(np.array([ar[:,1],ar[:,cols[i]]]))
    fits.HDUList([hdu]).writeto(name+'_'+perf+'.fits')

