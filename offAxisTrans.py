# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 09:43:17 2016

@author: cdelacroix
"""

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

lamD = fits.open(os.path.normpath(os.path.expandvars(path_lamD)))[0].data;
PSFs = fits.open(os.path.normpath(os.path.expandvars(path_PSFs)))[0].data;
offsets = [int(round(lamD[x]/lamDpp)) for x in range(lamD.size)];

# checking the background level and finding the centers
# pos 8 = 4.5 lambda/D = 143 pixels
# pos 14 = 9.0 lambda/D = 158 pixels
# 4.5 lambda/D = 15 pixels
bckg = [];
centers = [];
if False:
    for i in range(lamD.size): 
        print i
        bckg.append(np.median(PSFs[i,:,:]));
        (A,cx,cy,F) = fitFun.fit_airy_2D(PSFs[i,:,:]);
        centers.append([cx,cy])
    plt.figure(1); plt.plot(lamD,bckg);
    plt.xlabel(r"$\lambda$/D");plt.ylabel("Background level");
    plt.savefig("OAresults/bckg.png", dpi=300, transparent=True);plt.clf();
    plt.figure(2);
    plt.plot(lamD,zip(*centers)[0], label="x-direction");
    plt.plot(lamD,zip(*centers)[1], label="y-direction");
    plt.plot([12,0,12],[128,128,168],'k:');
    plt.xlabel(r"$\lambda$/D");
    plt.ylabel("2D-Airy-fitted center positions (pixels)");
    plt.legend(loc='upper left', fancybox=True, framealpha=0.3);
    plt.savefig("OAresults/centers.png", dpi=300, transparent=True);plt.clf();

# drawing all PSFs
myAxes = [-rim*lamDpp,rim*lamDpp,-rim*lamDpp,rim*lamDpp];
if False:
    for i in range(lamD.size): 
        psf = PSFs[i,cim-rim:cim+rim,cim-rim:cim+rim];
        plt.figure(3);plt.xlabel(r"$\lambda$/D");plt.ylabel(r"$\lambda$/D");
        plt.imshow(psf,origin='lower',extent=myAxes);plt.colorbar();
        plt.savefig("OAresults/_psf"+str(i)+".png", dpi=300, transparent=True);plt.clf();

# drawing sequence
sequence = np.array([0,2,3,4,5,10]);
fsz = 11; #fonsize
rim2 = int(2*FWHM);
vmin = 0; vmax = 1;
vmin = -12;vmax = 0;  # for log10
lim = (2*rim2+1)/2.*lamDpp;
if True:
    f = plt.figure(4);
    f.subplots_adjust(wspace=0.,right=.8);
    cbar_ax = f.add_axes([.81, 0.43, 0.01, 0.17]);
    for k,i in enumerate(sequence): 
        psf = PSFs[i,cim-rim2:cim+rim2+1,cim+offsets[i]-rim2:cim+offsets[i]+rim2+1]/np.max(PSFs);
        psf = np.log10(psf); # for log10
        myAxes = [-lim+lamD[i],lim+lamD[i],-lim,lim];
        ax = plt.subplot2grid((1,sequence.size),(0,k));
        im = ax.imshow(psf,origin='lower',extent=myAxes,vmin=vmin,vmax=vmax);
        ax.tick_params(labelsize=fsz,direction='out',length=3,right='off',top='off');
        ax.tick_params(left='off');
        ax.set_xticks([lamD[i]]);
        ax.set_yticks([-1,0,1]);  # is bugged, should appear on the figure...
        ax.add_artist(plt.Circle((lamD[i],0),0.5,color='w',fill=False));
        if k == 0:
            ax.tick_params(left='on');
            ax.set_ylabel(r"$\lambda$/D",fontsize=fsz)
            cb = f.colorbar(im, cax=cbar_ax);
            cb.ax.tick_params(labelsize=fsz);
            cb.set_ticks([vmin,(vmax+vmin)/2.,vmax]);
        if k == round(sequence.size/2):
            ax.set_xlabel(r"$\lambda$/D",fontsize=fsz);
    plt.text(5,2.3,"log10")
    plt.text(-5.5,2.3,"circles = FWHM")
    plt.setp([a.get_yticklabels() for a in f.axes[1:]], visible=False);   
    plt.savefig("OAresults/sequence.png", dpi=300, transparent=True);plt.clf();

# Radial profile
ioff = 10; # position of off-axis psf
xrim = [x*lamDpp for x in range(rim)];
if True:
    plt.figure(5);
    P_ofa = fitFun.get_radial_profile(PSFs[ioff],(cim+offsets[ioff],cim),nbin);
    P_ona = fitFun.get_radial_profile(PSFs[0],(cim+offsets[0],cim),nbin);
    plt.plot(xrim,P_ofa[:36]/np.max(P_ofa),label="Off-axis PSF");
    plt.plot(xrim,P_ona[:36]/np.max(P_ofa),label="On-axis PSF");
    plt.xlabel(r"$\lambda$/D");plt.ylabel("Radial profile");
    plt.yscale('log');plt.xlim(xrim[0],xrim[-1]);
    plt.legend(fancybox=True, framealpha=0.3);
    plt.savefig("OAresults/rad_prof.png", dpi=300, transparent=True);plt.clf();

# Circular aperture photometry
# radius of the integrated zone = FWHM/2 = 1/2 lambda/D
OAT = [];
circ = [];
ring = [];
IWA = 3.0;
OWA = 10.3;
if False:
    for i in range(lamD.size): 
        mask_circ = fitFun.get_r_dist(2*cim,2*cim,cim+offsets[i],cim) <= FWHM/2.;
        mask_ring = abs(fitFun.get_r_dist(2*cim,2*cim,cim,cim) - offsets[i]) <= FWHM/2.;
        OAT.append(np.sum(PSFs[i,:,:][mask_circ]));
        circ.append(np.mean(PSFs[i,:,:][mask_circ]));
        ring.append(np.mean(PSFs[0,:,:][mask_ring]))
    #off-axis transmission curve
    plt.figure(6);
    plt.plot(lamD,OAT);
    plt.ylim(1e-3,1e-1);plt.xlim(lamD[0],lamD[-1]);plt.yscale('log');
    plt.xlabel(r"$\lambda$/D");plt.ylabel("Off-axis transmission");
    plt.savefig("OAresults/OAT.png", dpi=300, transparent=True);#plt.clf();
    #normalized off-axis transmission curve
    plt.figure(7);
    plt.plot(lamD,circ/np.max(circ));
    plt.plot([lamD[0],lamD[-1]],[.5,.5],'k:');
    plt.plot([IWA,IWA],[.5,0],'k:');
    plt.plot([OWA,OWA],[.5,0],'k:');
    plt.text(3.5,.25," IWA = "+str(IWA)+" $\lambda$/D \nOWA = "+str(OWA)+" $\lambda$/D")
    plt.ylim(0,1);plt.xlim(lamD[0],lamD[-1]);
    plt.xlabel(r"$\lambda$/D");plt.ylabel("Normalized off-axis transmission");
    plt.savefig("OAresults/OATnorm.png", dpi=300, transparent=True);plt.clf();
    #contrast curve
    plt.figure(8);
    contrast = [ring[i]/circ[i] for i in range(lamD.size)];
    plt.plot(lamD,contrast);
    plt.yscale('log');plt.xlim(lamD[0],lamD[-1]);
    plt.xlabel(r"$\lambda$/D");plt.ylabel("Contrast");
    plt.savefig("OAresults/contrast.png", dpi=300, transparent=True);plt.clf();

    