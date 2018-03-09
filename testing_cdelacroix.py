# RUN ONE ENSEMBLE WITH SELECTED JSON SCRIPT
############################################

import os.path, EXOSIMS, EXOSIMS.MissionSim, pickle

# define the run_one function -- MUST specify the savefolder to save each DRM
def run_one(): 
    # wrap the run_sim in a try/except loop
    nbmax = 10
    for attempt in range(nbmax):
        try:
            # run one survey simulation
            SS.run_sim()
            res = SS.DRM[:]
        except:
            # if anything goes wrong, reset simulation
            sim.reset_sim()
        else:
            break
    else:
        raise ValueError("Unsuccessful run_sim after %s reset_sim attempts"%nbmax)
    # reset simulation at the end of each simulation
    SS.reset_sim(genNewPlanets=True, rewindPlanets=True)
    # need to pass savefolder to run_one
    import os.path
    filename = 'WFIRST_Nemati_KeplerLike.json'
    savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
    savefolder = os.path.join(savepath, filename[:-5]) #removes '.json'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    # store DRM of each survey
    import pickle
    from random import randint
    pklname = ''.join(["%s" % randint(0, 9) for num in range(9)]) + '.pkl'
    pklpath = os.path.join(savefolder, pklname)
    with open(pklpath, 'wb') as f:
        pickle.dump(res, f)
    # store OutSpec of each survey
    spcname = pklname[:-4] + '.spc'
    spcpath = os.path.join(savefolder, spcname)
    SS.genOutSpec(spcpath)
    return res


# number of survey simulations to run in the ensemble (e.g. 100 simulations)
nbsim = 2

# define scriptfile and savefolder
filename = 'WFIRST_Nemati_KeplerLike.json'
scriptpath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts')
scriptfile = os.path.join(scriptpath, filename)
savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
savefolder = os.path.join(savepath, filename[:-5])
if not os.path.exists(savefolder):
    os.makedirs(savefolder)

# run one Survey Ensemble -- ALSO SAVES the ensemble results
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
ens = sim.run_ensemble(nbsim, run_one=run_one)
enspath = os.path.join(savefolder, 'ensemble%s.pkl'%('_nbsim='+ str(int(nbsim))))
with open(enspath, 'wb') as f:
    pickle.dump(ens, f)


# CREATE A MAP OF SURVEY ENSEMBLES FOR MULTIPLE DMAG AND WA VALUES
#################################################################

import numpy as np
import os.path, EXOSIMS, EXOSIMS.MissionSim, pickle, time

# run_one does not save each DRM
def run_one(): 
    # wrap the run_sim in a try/except loop
    nbmax = 10
    for attempt in range(nbmax):
        try:
            # run one survey simulation
            SS.run_sim()
            res = SS.DRM[:]
        except:
            # if anything goes wrong, reset simulation
            sim.reset_sim()
        else:
            break
    else:
        raise ValueError("Unsuccessful run_sim after %s reset_sim attempts"%nbmax)
    # reset simulation at the end of each simulation
    SS.reset_sim(genNewPlanets=True, rewindPlanets=True)
    return res


# number of survey simulations to run in the ensemble (e.g. 100 simulations)
nbsim = 2

# define scriptfile and savefolder (with nbsim suffix)
filename = 'WFIRST_Nemati_KeplerLike.json'
scriptpath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts')
scriptfile = os.path.join(scriptpath, filename)
savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
savefolder = os.path.join(savepath, filename[:-5] + '_nbsim='+ str(int(nbsim)))
if not os.path.exists(savefolder):
    os.makedirs(savefolder)

# select MAP range values
varx = 'WAint'
varxvals = np.arange(.15,.43,.03)
#varxvals = np.arange(.15,.4,.04)
vary = 'dMagint'
varyvals = np.arange(20.,23.5,.5)
#varyvals = np.arange(20.5,23.5,.5)

nvals = len(varxvals)*len(varyvals)
t0 = time.time()
print '\nSimulation start: ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t0))
nfirst = 1  # change this value if the simulation crashed, to continue where it stopped
for i,x in enumerate(varxvals):
    for j,y in enumerate(varyvals):
        nval = i*len(varyvals)+j+1
        if nval >= nfirst:
            print '\nMap element #%s/%s' % (nval, nvals)
            print '%s = %s, %s = %s\n' % (varx, x, vary, y)
            reload(EXOSIMS)
            reload(EXOSIMS.MissionSim)
            sim = EXOSIMS.MissionSim.MissionSim(scriptfile, WAint=x, dMagint=y)
            ens = sim.run_ensemble(nbsim, run_one=run_one)
            enspath = os.path.join(savefolder, varx + '=' + str(x) + '_' + vary + '=' + str(y) + '.pkl')
            with open(enspath, 'wb') as f:
                pickle.dump(ens, f)
            tend = t0 + (time.time()-t0) / (nval-nfirst+1) * (nvals-nfirst+1)
            print '\nSimulation finish: ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(tend))


#########
#########


from EXOSIMS.util.analysis_tools import *
import EXOSIMS,EXOSIMS.MissionSim
import os.path,json
import numpy as np
%pylab

# LOAD MAP ENSEMBLES
nbsim = 30
varxvals = np.arange(.15,.43,.03)
#varxvals = np.arange(.15,.4,.04)
vary = 'dMagint'
varyvals = np.arange(20.,23.5,.5)
#varyvals = np.arange(20.5,23.5,.5)
savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
#filename = 'KnownRV_100d_0.4jit.json'
filename = 'KeplerLike_100d_1.6jit.json'
#filename = 'KeplerLike_1yr_full_MAP.json'
savefolder = savepath + '/' + filename[:-5] + '_nbsim='+ str(int(nbsim))
det_mean = np.zeros((len(varxvals),len(varyvals)))
det_std = np.zeros((len(varxvals),len(varyvals)))
status = 1
for i,x in enumerate(varxvals):
    for j,y in enumerate(varyvals):
        pth = savefolder + '/' + varx + '=' + str(x) + '_' + vary + '=' + str(y)
        ens = load_obj(pth)
        subfolder = savefolder + '/' + str(j)
        det = det_ensemble2(ens, status)
        det_mean[i,j] = mean(det)
        det_std[i,j] = std(det)


myAxes = [varxvals[0],varxvals[-1],varyvals[0],varyvals[-1]]
fig, ax = plt.subplots(figsize=(6,6))
#im = ax.imshow(det_mean,origin='lower',extent=myAxes)
im = ax.imshow(det_std,origin='lower',extent=myAxes)
ax.set_aspect(.1)
xlabel('WA (arcsec)')
ylabel('dMag')
fig.colorbar(im)
#title(filename[:-5] + ' -- MEAN')
#savefig(filename[:-5]+"_MEAN.png", dpi=300, transparent=True)
title(filename[:-5] + ' -- STD')
savefig(filename[:-5]+"_STD.png", dpi=300, transparent=True)



# RUN ONE ENSEMBLE OF SIMULATIONS
# -------------------------------
nbsim = 100
savefolder = savepath + filename[:-5]
savefolder = os.path.normpath(os.path.expandvars(savefolder))
# Check if save folder exists and create if necessary
if not os.path.exists(savefolder):
    os.makedirs(savefolder)
#genOutSpec_ensemble(scriptfile, savefolder, nbsim=nbsim)


# Load detection results from ensemble of simulations
unq = det_ensemble(savefolder)
det = det_ensemble(savefolder, status=1)
md = det_ensemble(savefolder, status=0)
out_iwa = det_ensemble(savefolder, status=-1)
out_owa = det_ensemble(savefolder, status=-2)
# Load characterization results from ensemble of simulations
full_spec = char_ensemble(savefolder, status=1)
part_spec = char_ensemble(savefolder, status=-1)
no_char = char_ensemble(savefolder, status=0)


# draw PDF curves
#draw_pdf(det, xlab='Total planet detections and MD')
draw_pdf(var1=det, var2=md, label1='detections', label2='missed detections', \
        xlab='Total number of detections',xmax=35)
savefig("det_pdf.png", dpi=300, transparent=True)
draw_pdf(var1=out_iwa, var2=out_owa, label1='below IWA', label2='beyond OWA', \
        xlab='Planets out of IWA-OWA range',xmax=55)
savefig("iwaowa_pdf.png", dpi=300, transparent=True)
draw_pdf(var1=full_spec, var2=part_spec, label1='Full specta', label2='Partial spectra', \
        xlab='Total number of characterizations',xmax=7)
savefig("char_pdf.png", dpi=300, transparent=True)



# RUN MULTIPLE ENSEMBLES OF SIMULATIONS
# -------------------------------------
varname = 'dMagint'
varvalues = np.arange(20,23.5,.5)
#varname = 'WAint'
#varvalues = np.arange(.15,.43,.03)
#varname = 'magEZ'
#varvalues = np.arange(20,24.5,.5)
nbsim = 30
savefolder = savepath + filename[:-5] + '_' + varname
savefolder = os.path.normpath(os.path.expandvars(savefolder))
for var in varvalues:
    subfolder = savefolder + '/' + str(var)
    # Check if save folder exists and create if necessary
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    genOutSpec_ensemble(scriptfile, subfolder, nbsim=nbsim)


# Load mean detection results of multiple ensembles
#mean_det = multi_det(savefolder, status=1)


draw_multi_ensemble(savefolder, varvalues, status=None, label1='mean', label2='std', xlab=varname+' [arcsec]', \
        ylab='Unique Detections', linewidth=2, framealpha=0.3, markersize=5, fontsize=18)
savefig("multi_unq.png", dpi=300, transparent=True)
draw_multi_ensemble(savefolder, varvalues, status=1, label1='mean', label2='std', xlab=varname+' [arcsec]', \
        ylab='Total Detections', linewidth=2, framealpha=0.3, markersize=5, fontsize=18)
savefig("multi_det.png", dpi=300, transparent=True)
draw_multi_ensemble(savefolder, varvalues, status=0, label1='mean', label2='std', xlab=varname+' [arcsec]', \
        ylab='Missed Detections', linewidth=2, framealpha=0.3, markersize=5, fontsize=18)
savefig("multi_md.png", dpi=300, transparent=True)
draw_multi_ensemble(savefolder, varvalues, status=-1, label1='mean', label2='std', xlab=varname+' [arcsec]', \
        ylab='Below IWA', linewidth=2, framealpha=0.3, markersize=5, fontsize=18)
savefig("multi_iwa.png", dpi=300, transparent=True)
draw_multi_ensemble(savefolder, varvalues, status=-2, label1='mean', label2='std', xlab=varname+' [arcsec]', \
        ylab='Beyond OWA', linewidth=2, framealpha=0.3, markersize=5, fontsize=18)
savefig("multi_owa.png", dpi=300, transparent=True)








