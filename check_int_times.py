import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u

filename = 'WFIRST_Nemati_KeplerLike.json'
scriptpath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts')
#scriptpath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/Scripts')
scriptfile = os.path.join(scriptpath, filename)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)

SS = sim.SurveySimulation
SU = sim.SimulatedUniverse
TL = sim.TargetList
OS = sim.OpticalSystem
Obs = sim.Observatory
ZL = sim.ZodiacalLight
TK = sim.TimeKeeping

# selection criteria
nStars = 50
Rp0 = 2*u.earthRad
sma0 = 2*u.AU

# find planets and stars of interest
err = np.abs(SU.Rp - Rp0).to('earthRad').value + np.abs(SU.a - sma0).to('AU').value
pInds = np.array([], dtype=int)
sInds = np.array([], dtype=int)
cnt = 0
while cnt < nStars:
    pInd = np.argmin(err)
    sInd = SU.plan2star[pInd]
    if sInd not in sInds:
        pInds = np.append(pInds, pInd)
        sInds = np.append(sInds, sInd)
        cnt += 1
    err[pInd] = np.inf
# sort indices
pInds.sort()
sInds.sort()

# update Simulated Universe
SU.revise_planets_list(pInds)
SU.revise_stars_list(sInds)
sim.reset_sim(genNewPlanets=False)
print SU.plan2star

# run simulation
sim.run_sim()

#######
#######

# comparing det_time with minimum integration times for detection (at SNR=5)

specs = json.loads(open(os.listdir('./')[0]).read())

bug = []
cnt = 0
#while not(bug):
folder = '../testing/results/wfirst_nemati_nom_6yr_full2/'
names = [x for x in os.listdir('./') if x[-4:] == '.pkl' and x not in '.DS_Store']
for name in names:
    with open(name, 'rb') as f:
        res = pickle.load(f)
#    sim.reset_sim()
#    sim.run_sim()
#    res = {}
#    res['DRM'] = SS.DRM
#    res['systems'] = SU.dump_systems()
    det_t = np.hstack([np.array([row['det_time'].value]*np.sum(row['det_status'] == 1)) for row in res['DRM']])*row['det_time'].unit
    det_SNR = np.hstack([row['det_SNR'][row['det_status'] == 1] for row in res['DRM']])
    det_fZ = np.hstack([np.array([row['det_fZ'].value]*np.sum(row['det_status'] == 1)) for row in res['DRM']])*row['det_fZ'].unit
    det_fEZ = np.hstack([row['det_params']['fEZ'][row['det_status'] == 1].value for row in res['DRM']])*row['det_params']['fEZ'].unit
    det_dMag = np.hstack([row['det_params']['dMag'][row['det_status'] == 1] for row in res['DRM']])
    det_WA = np.hstack([row['det_params']['WA'][row['det_status'] == 1].value for row in res['DRM']])*row['det_params']['WA'].unit
    mode = OS.observingModes[0]
    pInds = np.hstack([row['plan_inds'][row['det_status'] == 1] for row in res['DRM']]).astype(int)
    sInds = res['systems']['plan2star'][pInds]
    det_t_min = OS.calc_intTime(TL, sInds, det_fZ, det_fEZ, det_dMag, det_WA, mode)
    
    bugi = np.any(det_t - det_t_min < 0)
    cnt += 1
    print '\n Survey #%s: bug is %s \n'%(cnt, bugi)
    bug.append(bugi)
