import os.path
import glob
import EXOSIMS, EXOSIMS.MissionSim
import cPickle as pickle
import astropy.units as u
import numpy as np
import matplotlib.pylab as plt

from miscpy.PlotFun.logHist import logHist



def gen_summary(run_name):
    basedir = 'CODE/Python/EXOSIMS/testing/results/'
    pklfiles = glob.glob(os.path.join(basedir,run_name,'*.pkl'))

    detected = []
    fullspectra = []
    partspectra = []
    Rps = []
    Mps = []
    tottime = []
    starinds = []
    smas = []

    for f in pklfiles:
        print f
        with open(f, 'rb') as g:
            res = pickle.load(g)

        dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
        detected.append(dets)
        fullspectra.append(np.hstack([np.array(row['plan_inds'])[np.array(row['char_status']) == 1] for row in res['DRM']]))
        partspectra.append(np.hstack([np.array(row['plan_inds'])[np.array(row['char_status']) == -1] for row in res['DRM']]))
        tottime.append(np.sum([row['det_time']+row['char_time'] for row in res['DRM']]))
        Rps.append(res['systems']['Rp'][np.unique(dets)].to('earthRad').value)
        smas.append(res['systems']['a'][np.unique(dets)].to('AU').value)
        Mps.append(res['systems']['Mp'][np.unique(dets)].to('earthMass').value)
        starinds.append(np.array([row['star_ind'] for row in res['DRM']]))

    return detected, fullspectra,partspectra,Rps, Mps,tottime,starinds,smas


Rpin = []
basedir = '/data/extmount/EXOSIMSres'
pklfiles = glob.glob(os.path.join(basedir,'wfirst_nemati_run8','*.pkl'))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
        Rpin.append(res['systems']['Rp'].to('earthRad').value)





d1,f1,p1,Rp1,Mp1,tt1 = gen_summary('wfirst_nemati_run1') #checkkeepout=false, dMagLim =22.5, bad
d2,f2,p2,Rp2,Mp2,tt2 = gen_summary('wfirst_nemati_run2') #checkkeeoput=true, bad
d3,f3,p3,Rp3,Mp3,tt3 = gen_summary('wfirst_nemati_run3') #dMagint = 23, bad

d4,f4,p4,Rp4,Mp4,tt4,sts4,smas4 = gen_summary('wfirst_nemati_run4') #dMagint = 22.5,intcutoff = 30
d5,f5,p5,Rp5,Mp5,tt5 = gen_summary('wfirst_nemati_run5') #dMagint = 23,intcutoff = 30
d6,f6,p6,Rp6,Mp6,tt6 = gen_summary('wfirst_nemati_run6') #dMagint = 22.5, intCutoff = 60
d7,f7,p7,Rp7,Mp7,tt7 = gen_summary('wfirst_nemati_run7') #intCutoff = 100
d8,f8,p8,Rp8,Mp8,tt8 = gen_summary('wfirst_nemati_run8') #intCutoff = 30, dMagint = 22
d9,f9,p9,Rp9,Mp9,tt9 = gen_summary('wfirst_nemati_run9') #intCutoff = 30, dMagint = 21

d10,f10,p10,Rp10,Mp10,tt10 = gen_summary('wfirst_nemati_run10') #intCutoff = 30, dMagint = 21




def dist_plot(res,uniq = True,fig=None,lstyle='--',legtext=None):
    rcounts = []
    for el in res:
        if uniq:
            rcounts.append(np.array([np.unique(r).size for r in el]))
        else:
            rcounts.append(np.array([len(r) for r in el]))

    bins = range(np.min(np.hstack(rcounts)),np.max(np.hstack(rcounts))+2)

    syms = 'osp^v<>h'
    if fig is None:
        plt.figure()
    else:
        plt.figure(fig)
        plt.gca().set_prop_cycle(None)

    if legtext is None:
        legtext = [None]*len(res)

    for j in range(len(res)):
        tmp = np.histogram(rcounts[j],bins=bins)[0].astype(float)
        plot(bins[:-1], tmp/tmp.max(), syms[np.mod(j,len(syms))]+lstyle,label=legtext[j])



dist_plot([d4,d8,d10],legtext = ['$\\Delta$mag=22.5','$\\Delta$mag=22.0','$\\Delta$mag=21.0'],lstyle='-') 
dist_plot([f4,f8,f10],fig=1) 



#-------------------------------------------
logHist(np.hstack(Rpin).value,fig=1,pdf=True,linestyle='--',color='k',label='Input')
logHist(np.hstack(Rp4).value,fig=1,pdf=True,noclear=True,label='$\\Delta$mag=22.5')
logHist(np.hstack(Rp10).value,fig=1,pdf=True,noclear=True,label='$\\Delta$mag=22')
plt.xlabel('R ($R_\\oplus$)')
plt.ylabel('$f_\\bar{R}(R)$')
plt.legend()





with open('./EXOSIMSres/data.pkl', 'wb') as g:
    pickle.dump({'run4':(d4,f4,p4,Rp4,Mp4,tt4),'run5':(d5,f5,p5,Rp5,Mp5,tt5),'run6':(d6,f6,p6,Rp6,Mp6,tt6),'run7':(d7,f7,p7,Rp7,Mp7,tt7),'run10':(d10,f10,p10,Rp10,Mp10,tt10),'run8':(d8,f8,p8,Rp8,Mp8,tt8),'Rpin':Rpin},g)


#=======================
aedges = np.arange(0,13,0.25)
Redges = np.arange(0,17,0.25)


h = np.zeros((aedges.size-1,Redges.size-1))
for R,a in zip(Rp4,smas4):
    h += np.histogram2d(a,R,bins=[aedges,Redges])[0]

acents = np.diff(aedges)/2.+aedges[:-1]
Rcents = np.diff(Redges)/2.+Redges[:-1]


hin = np.histogram2d(res['systems']['a'].value,(res['systems']['Rp']/u.R_earth).decompose().value,bins=[30,20])


# global input
nx = 100
ny = 30
basedir = 'CODE/Python/EXOSIMS/testing/results/'
#run_name = 'wfirst_nemati_run4'
run_name = 'wfirst_nemativ2_run1'
pklfiles = glob.glob(os.path.join(basedir, run_name, '*.pkl'))


# radius vs sma
###############
xedges = np.logspace(-1, 2, nx+1)
yedges = np.logspace(0, np.log10(22.6), ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    hin += np.histogram2d(res['systems']['a'].to('AU').value,
            res['systems']['Rp'].to('earthRad').value,
            bins=[xedges, yedges])[0]
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    h += np.histogram2d(res['systems']['a'][np.unique(dets)].to('AU').value,
            res['systems']['Rp'][np.unique(dets)].to('earthRad').value,
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
#im = plt.pcolor(xedges, yedges, np.log10(h.T))
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
im.cmap.set_under('w')
#ax.set_aspect(2)
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$R \, / \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)


# eccen vs sma
##############
xedges = np.logspace(-1, 2, nx+1)
yedges = np.linspace(0, 1, ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    # select planets with radius = [1.5-2.5] earthRad
    mask = np.abs(res['systems']['Rp'][np.unique(dets)].to('earthRad').value - 2) <= 0.5
    h += np.histogram2d(res['systems']['a'][np.unique(dets)][mask].to('AU').value,
            res['systems']['e'][np.unique(dets)][mask],
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
im.cmap.set_under('w')
#ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$e$')
plt.title('$R = 1.5 - 2.5 \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)

# tint vs sma
##############
xedges = np.logspace(-1, 2, nx+1)
yedges = np.logspace(-2, 2, ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
h_R = np.zeros((nx, ny))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    tdets = np.hstack([np.array([row['det_time']]*np.sum(np.array(row['det_status']) == 1)) for row in res['DRM']])
    h += np.histogram2d(res['systems']['a'][dets].to('AU').value,
            tdets,
            bins=[xedges, yedges])[0]
    # select planets with radius = [1.5-2.5] earthRad
    mask = np.abs(res['systems']['Rp'][dets].to('earthRad').value - 2) <= 0.5
    h_R += np.histogram2d(res['systems']['a'][dets][mask].to('AU').value,
            tdets[mask],
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
#im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h_R.T), cmap=cmaps.viridis)
#im.cmap.set_under('w')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$detection \, t_{int} \, (day)$')
#plt.title('All radii')
plt.title('$R = 1.5 - 2.5 \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)


# WA vs sma
###############
xedges = np.logspace(-1, 2, nx+1)
yedges = np.logspace(2, 3.5, ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
h_R = np.zeros((nx, ny))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    WAs = np.hstack([np.array(row['det_WA'])[np.array(row['det_status']) == 1] for row in res['DRM'] if np.any(row['plan_inds'])])
    h += np.histogram2d(res['systems']['a'][dets].to('AU').value,
            WAs,
            bins=[xedges, yedges])[0]
    # select planets with radius = [1.5-2.5] earthRad
    mask = np.abs(res['systems']['Rp'][dets].to('earthRad').value - 2) <= 0.5
    h_R += np.histogram2d(res['systems']['a'][dets][mask].to('AU').value,
            WAs[mask],
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
#im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h_R.T), cmap=cmaps.viridis)
#im.cmap.set_under('w')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$WA \, (mas)$')
#plt.title('All radii')
plt.title('$R = 1.5 - 2.5 \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)


# fEZ vs sma
###############
xedges = np.logspace(-1, 2, nx+1)
yedges = np.logspace(-11, -6, ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
h_R = np.zeros((nx, ny))
for f in pklfiles:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    fEZs = np.hstack([np.array(row['det_fEZ'])[np.array(row['det_status']) == 1] for row in res['DRM'] if np.any(row['plan_inds'])])
    h += np.histogram2d(res['systems']['a'][dets].to('AU').value,
            fEZs,
            bins=[xedges, yedges])[0]
    # select planets with radius = [1.5-2.5] earthRad
    mask = np.abs(res['systems']['Rp'][dets].to('earthRad').value - 2) <= 0.5
    h_R += np.histogram2d(res['systems']['a'][dets][mask].to('AU').value,
            fEZs[mask],
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
#im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h_R.T), cmap=cmaps.viridis)
#im.cmap.set_under('w')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$fEZ \, (1/arcsec2)$')
plt.title('All radii')
#plt.title('$R = 1.5 - 2.5 \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)



# check int times
#################
xedges = np.logspace(-1, 2, nx+1)
yedges = np.logspace(-2, 2, ny+1)
hin = np.zeros((nx, ny))
h = np.zeros((nx, ny))
h_R = np.zeros((nx, ny))
for f in pklfiles[:2]:
    print f
    with open(f, 'rb') as g:
        res = pickle.load(g)
    tdets = np.hstack([np.array([row['det_time']]*np.sum(np.array(row['det_status']) == 1)) for row in res['DRM']])*u.day
    dets = np.hstack([np.array(row['plan_inds'])[np.array(row['det_status']) == 1] for row in res['DRM']]).astype(int)
    sInds = res['systems']['plan2star'][dets]
    Rp = res['systems']['Rp'][dets]

    mode = OS.observingModes[0]
    fZ = np.hstack([np.array([row['det_fZ']]*np.sum(np.array(row['det_status']) == 1)) for row in res['DRM']])/u.arcsec**2
    fEZ = np.hstack([np.array(row['det_fEZ'])[np.array(row['det_status']) == 1] for row in res['DRM'] if np.any(row['plan_inds'])])/u.arcsec**2
    dMag = np.hstack([np.array(row['det_dMag'])[np.array(row['det_status']) == 1] for row in res['DRM'] if np.any(row['plan_inds'])])
    WA = np.hstack([np.array(row['det_WA'])[np.array(row['det_status']) == 1] for row in res['DRM'] if np.any(row['plan_inds'])])*u.mas
    tdets_min = OS.calc_intTime(TL, sInds, fZ, fEZ, dMag, WA, mode)
    C_p, C_b, C_sp = OS.Cp_Cb_Csp(TL, sInds, fZ, fEZ, dMag, WA, mode)
    # calculate signal and noise levels (based on Nemati14 formula)
    Signal = (C_p*tdets_min).decompose().value
    Noise = np.sqrt((C_b*tdets_min + (C_sp*tdets_min)**2).decompose().value)
    SNR = Signal/Noise

    h += np.histogram2d(res['systems']['a'][dets].to('AU').value,
            tdets_min.to('day').value,
            bins=[xedges, yedges])[0]
    # select planets with radius = [1.5-2.5] earthRad
    mask = np.abs(res['systems']['Rp'][dets].to('earthRad').value - 2) <= 0.5
    h_R += np.histogram2d(res['systems']['a'][dets][mask].to('AU').value,
            tdets_min[mask].to('day').value,
            bins=[xedges, yedges])[0]

import colormaps as cmaps
plt.figure()
fig, ax = plt.subplots(figsize=(6,6))
im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h.T), cmap=cmaps.viridis)
#im = plt.contourf(xedges[:-1], yedges[:-1], np.log10(h_R.T), cmap=cmaps.viridis)
#im.cmap.set_under('w')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('$a \, (AU)$')
plt.ylabel('$t_{int}$')
plt.title('All radii')
#plt.title('$R = 1.5 - 2.5 \, R_\\oplus$')
fig.colorbar(im)
plt.savefig("grid.png", dpi=300, transparent=True)
