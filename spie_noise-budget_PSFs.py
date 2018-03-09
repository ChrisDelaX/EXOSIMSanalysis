%pylab

import os.path
import numpy as np
import matplotlib.pyplot as plt  
from astropy.io import fits 
import fitFun #lib with fitting functions (Gauss, Airy,...) and radial profile

# Input Data
# ----------
instru = 'WFIRST';
path_lamD = "$HOME/INSTRUMENTS/WFIRST/offaxis_psfs/psf_offset_lamDivD.fits";
path_PSFs = "$HOME/INSTRUMENTS/WFIRST/offaxis_psfs/hlc20140623-139_offaxis_psfs.fits";
minprofiles = 1e-4; # minimal radial value, to allow log scale
nbin = 1.; # radial profile binning parameter
lamDpp = 0.3; #lambda/D per pixel
cim = int(128); #image center (pix)
rim = int(36); #image radius (pix)
FWHM = 1/lamDpp; # Full Width at Half Maximum = lambda/D
myAxes = [-rim*lamDpp,rim*lamDpp,-rim*lamDpp,rim*lamDpp];
rim2 = int(2*FWHM);

lamD = fits.open(os.path.normpath(os.path.expandvars(path_lamD)))[0].data;
PSFs = fits.open(os.path.normpath(os.path.expandvars(path_PSFs)))[0].data;
offsets = [int(round(lamD[x]/lamDpp)) for x in range(lamD.size)];

sMag = 5
dMag = 22.5
#vmin=-9;vmax=-11.5;
vmin=-8.1;vmax=-8.4;
sInd = 0

pPSF1 = 1/np.max(PSFs)*PSFs[2,cim-rim:cim+rim,cim-rim:cim+rim]*10.**(-0.4*dMag) # 2.1 lamoverd
pPSF2 = 1/np.max(PSFs)*rot90(PSFs[4,cim-rim:cim+rim,cim-rim:cim+rim],1)*10.**(-0.4*dMag) # 2.7 lamoverd
pPSF3 = 1/np.max(PSFs)*rot90(PSFs[7,cim-rim:cim+rim,cim-rim:cim+rim],2)*10.**(-0.4*dMag) # 3.9 lamoverd
pPSF4 = 1/np.max(PSFs)*rot90(PSFs[11,cim-rim:cim+rim,cim-rim:cim+rim],3)*10.**(-0.4*dMag) # 6 lamoverd
sPSF = 1/np.max(PSFs)*PSFs[sInd,cim-rim:cim+rim,cim+offsets[sInd]-rim:cim+offsets[sInd]+rim]
zodi = np.ones(shape(pPSF))*np.max(pPSF)/5
cam = np.ones(shape(pPSF))*np.max(pPSF)*5

#vmin=0;vmax=-2.6;
#sInd = 6
#sPSF = 1/np.max(PSFs[sInd,:,:])*PSFs[sInd,cim-rim:cim+rim,cim+offsets[sInd]-rim:cim+offsets[sInd]+rim]


nP = 1#1.
nZ = 1#1.
nC = 1#1.
PP = 1./30
close()
plot(0.1,0.1,'w*',markersize=16)
imshow(log10(nP*(pPSF1+pPSF2+pPSF3+pPSF4)+PP*sPSF+nZ*zodi+nC*cam),origin='lower',extent=myAxes,vmin=vmin,vmax=vmax)
colorbar();

savefig("psf.png", dpi=300, transparent=True);
