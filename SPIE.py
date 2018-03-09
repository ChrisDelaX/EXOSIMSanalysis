import os.path,EXOSIMS,EXOSIMS.MissionSim

import os, inspect, copy, numbers, json, pandas, h5py 
import numpy as np
import cPickle as pickle
import astropy.io.fits as fits
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d, interp2d, griddata
from EXOSIMS.util.deltaMag import deltaMag
from EXOSIMS.util.get_module import get_module
from EXOSIMS.util.eccanom import eccanom

folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/exomissionsim/MyFiles/Scripts'))
for i in range(1):
    print 'Loop index', i
    reload(EXOSIMS)
    reload(EXOSIMS.MissionSim)
#    filename = 'template_WFIRST_KnownRV.json'
#    filename = 'template_WFIRST_EarthTwinHabZone.json'
#    filename = 'HabX_4m_HZ_20160621_1114.json'
    filename = 'HabX_4m_HZ_20160719nopara.json'
    scriptfile = os.path.join(folder,filename)
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
    
    Obs = sim.Observatory
    OS = sim.OpticalSystem
    TL = sim.TargetList
    TK = sim.TimeKeeping
    SU = sim.SimulatedUniverse
    PPro = sim.PostProcessing
    PPop = sim.PlanetPopulation
    PPMod = sim.PlanetPhysicalModel
    BS = sim.BackgroundSources
    ZL = sim.ZodiacalLight
#    specs = json.loads(open(scriptfile).read())
#    SC = get_module(specs['modules']['StarCatalog'],'StarCatalog')(**specs)
    SS = sim.SurveySimulation
    SS.run_sim()
    SS._outspec['DRM'] = SS.DRM
#    sim.genOutSpec('Results/Nemati_' + str(sim._outspec['seed']))


#import cProfile as profile
#profile.run('sim.SurveySimulation.run_sim()','output.qcg')


DRM = sim.SurveySimulation.DRM
#number of observations
nobs = len(DRM)
#detection status
ds = [DRM[i]['det_status'] for i in range(nobs)]
det = np.array([int(ds[i] >0) for i in range(nobs)])
#number of detections
ndet = np.sum(det)
#characterization status
cs = [DRM[i].get('char_1_time',np.nan) for i in range(nobs)]
char = np.array([int(cs[i] >0) for i in range(nobs)])
#number of characterizations
nchar = np.sum(char)
#planets detected for each observation
pInds = [DRM[i]['plan_inds'] for i in range(nobs)]
#working angles of these planets
pWA = np.array([DRM[i].get('det_WA',np.nan) for i in range(nobs)])*u.mas
#creating an array with unique planet detections
vp = []
wp = []
for v in pInds:
    if v:
        vp.append(v) # array-of-arrays of detected planets
    for w in v:
        wp.append(w) # array of detected planets
pUnique = np.unique(wp) 


np.where(det)[0]
np.where(char)[0]
print TL.nStars, SU.nPlans
print nobs, pInds, ndet, nchar
