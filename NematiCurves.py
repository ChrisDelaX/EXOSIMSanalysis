#%pylab
#from matplotlib.pyplot import *
import numpy as np


r_CG_i = np.array([10**-4.8,10**-4.5,10**-4.1,10**-3.8,10**-3.8,10**-3.6,\
        10**-3.3,10**-3.3,10**-3.1,10**-3.1,10**-2.8,10**-2.8,10**-2.8,\
        10**-2.8,10**-2.6,10**-2.6,10**-2.3,10**-2.3,10**-2.3,10**-2.3,\
        10**-2.1,10**-1.9,10**-1.9,10**-1.4,10**-1.1,10**-0.2]);
#shuffle(r_CG_i)
r_CG_i = r_CG_i[6:22]  # [:16] [15::-1]

Q_pl_i = np.array([10**-9.0,10**-8.9,10**-8.9,10**-8.9,10**-8.7,10**-8.6,\
        10**-8.6,10**-8.5,10**-8.3,10**-8.2,10**-8.0,10**-8.0,10**-8.0,\
        10**-7.9,10**-7.8,10**-7.8]);

# Q_pl_med = 3e-9; #contrast
# r_CG_med = 3e-3; #(e-/s)

D = 2.4;
s = np.pi/4.;
A = s*D**2;#3.5#
Apixel = np.pi*(0.7*660e-9/D*180*3600/np.pi)**2;

tau = 0.01#0.001742;#
QE = 0.88;  # QE at 660, 770, 890 nm (88, 68, 28%)  # 0.88 / 0.68 / 0.52
Vmag = 5;
F0 = 3631*1.51e7;
BW = 0.20;#.18#
r_psf = BW*F0*10**(-Vmag/2.5)*QE*tau*A;
# contrast before post processing
f_PP = 0.1;
Q_PP = np.median(r_CG_i)/r_psf;
Q = Q_PP/f_PP;
print Q         # 4.7e-9 / 6.1e-9 / 7.9e-9

# general parameters
texp = 1000;
mpix = 4;#22.3#
Nspec = 70*BW#14;
f_SR = 1./Nspec;#1.#

# different cameras
cam = np.array(['CCD','AM-EMCCD','PC-EMCCD','sCMOS']);
ENF = np.array([1., np.sqrt(2), 1., 2.]);
q_cic = np.array([1e-3, 1e-3, 1e-3, 0.]);
read_noise = np.array([3., 16./500, 16./500, 1.]);

r_pl = f_SR*r_psf * Q_pl_i;  
r_CG = f_SR*r_CG_i;                # for Bijan's curves
#r_CG = f_SR*r_psf *f_PP*6e-9
r_zodi = f_SR*r_psf * 10**(-23.54/2.5)*8*Apixel;
r_dark = 5e-4*mpix;
r_cic = q_cic*mpix/texp;
r_read = (read_noise)**2*mpix/texp;

# calculate integration time
SNR = 5
for i in range(3):
    r_noise = ENF[i]**2*(r_pl + r_CG + r_zodi + r_dark + r_cic[i]) + r_read[i];
    t = (SNR**2 * r_noise) / (r_pl**2 - SNR**2 * r_CG**2);
    t = t[t>0]/(24*3600);
    t.sort();
    x = range(len(t)+1)[1:];
    y = np.empty(len(t))
    for j in range(len(t)):
        if j == 0:
            y[j]=t[j];
        else:
            y[j]=t[j]+y[j-1];
    plot(x,y,'-*',label=cam[i]);

grid('on');xlabel('Number of candidates');ylabel('Cumulative time (days)')
plot([1,15],[30,30],'k--')
yscale('log');xlim(1,15);ylim(6e-2,6e1)
legend(fancybox=True, framealpha=0.3, loc='lower right');


