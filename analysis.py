# python exosims-run EXOSIMS/Scripts/sampleScript_coron.json
# ipython EXOSIMS/e2eTests.py
# rsync -av INSTRUMENTS/WFIRST/ cdelacroix@atuin.coecis.cornell.edu:INSTRUMENTS/WFIRST/
# rsync -av INSTRUMENTS/HabEx/ cdelacroix@atuin.coecis.cornell.edu:INSTRUMENTS/HabEx/
# rsync -av CODE/Python/EXOSIMS/ cdelacroix@atuin.coecis.cornell.edu:CODE/Python/EXOSIMS/
# rsync -av cdelacroix@atuin.coecis.cornell.edu:CODE/Python/EXOSIMS/ CODE/Python/EXOSIMS/
# rsync -av cdelacroix@atuin.coecis.cornell.edu:../../data/extmount/EXOSIMSres/wfirst_nemativ2_run1 CODE/Python/EXOSIMS/testing/results/wfirst_nemativ2_run1
# rsync -av cdelacroix@atuin.coecis.cornell.edu:../../data/extmount/EXOSIMSres/wfirst_nemati_nom_6yr_full2 CODE/Python/EXOSIMS/testing/results/wfirst_nemati_nom_6yr_full2

# ssh cdelacroix@atuin.coecis.cornell.edu
# cd CODE/Python/EXOSIMS/gitexosims/
# ipcluster start -n 32
# python partest.py

# import os.path,EXOSIMS,EXOSIMS.MissionSim
# scriptfile = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts/sampleScript_coron.json')
# sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
# sim.run_sim()
# 
# DRM = sim.SurveySimulation.DRM
# sInd = [DRM[x]['star_ind'] for x in range(len(DRM))]
# ra = sim.TargetList.coords.ra[sInd]
# dec = sim.TargetList.coords.dec[sInd]
# %pylab
# plot(ra,dec)
# axis([0, 360, -90, 90])


import os.path, EXOSIMS, EXOSIMS.MissionSim
import os, inspect, copy, numbers, json, pandas, logging, pickle, subprocess, itertools, timeit
import h5py
import numpy as np
import cPickle as pickle
import astropy.io.fits as fits
import astropy.units as u
import astropy.constants as const
import scipy.stats as st
from scipy.interpolate import interp1d, interp2d, griddata
from astropy.time import Time
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
from EXOSIMS.util.deltaMag import deltaMag
from EXOSIMS.util.get_module import get_module
from EXOSIMS.util.eccanom import eccanom
from EXOSIMS.util import statsFun 
Logger = logging.getLogger(__name__)

#folder = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts')
folder = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/Scripts')
for i in range(1):
    print 'Loop index', i
    reload(EXOSIMS);reload(EXOSIMS.MissionSim);
#    filename = 'outspec.json'
#    filename = 'WFIRST_Nemati_KnownRV.json'
#    filename = 'WFIRST_Nemati_KeplerLike.json'
#    filename = 'HabEx_4m_SS_20170411.json'
    filename = 'Rhonda/HabEx_4m_SAG13_cluster_20170613_0805.json'
#    filename = 'Rhonda/HabEx_4m_TS_20170727.json'
#    filename = 'Rhonda/HabEx_4m_SS_232_cluster_20170804.json'
#    filename = 'Rhonda/HabEx_4m_SAG13_20170610.json'
#    filename = 'sampleScript_coron.json'
#    filename = 'parscript.json'
#    filename = 'TestScripts/01_all_defaults.json'
#    filename = 'TestScripts/02_KnownRV_FAP=1_WFIRSTObs_staticEphem.json'
#    filename = 'TestScripts/03_EarthTwin_Coronagraph_GarrettComp.json'
#    filename = 'TestScripts/04_KeplerLike_Occulter_linearJScheduler.json'
#    filename = 'KeplerLike_cbyt_100d_0.4jit.json'
#    filename = 'KeplerLike_1yr_full_MAP.json'
#    filename = 'KnownRV_1yr_full_MAP.json'
#    filename = 'walker-dula.json'
#    filename = 'HabExAYO_4m_local_20170411.json'
#    filename = 'AYO_8m_20161020.json'
#    filename = 'sS_AYO.json'
#    filename = 'mike_bottom.json'
#    filename = 'rahul_template_WFIRST_KnownRV.json'
#    filename = 'template_rpateltest_KnownRV_2years.json'
#    filename = 'HabX_4m_400nm_20161107.json'
#    filename = 'HabX_8m_400nm_20160927_test.json'
#    filename = 'SAG13_4m_400nm_20161202.json'
#    filename = 'WFIRST_KnownRV_forComparison.json'
    scriptfile = os.path.join(folder, filename)
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile, keepStarCatalog=True)
    Comp = sim.Completeness
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
    SS = sim.SurveySimulation


##

    sim.run_sim()



#######
#######


# class Object(object):
#     def __init__(self,  **specs):
#         for att in specs:
#             setattr(self,att,lambda x=0: specs[att])
# 
#     def test(self):
#         print 10
# 
# self = Object(test=5)
# 
# print self.test()




#    specs = json.loads(open(scriptfile).read())
#    SC = get_module(specs['modules']['StarCatalog'],'StarCatalog')(**specs)
#    sim.genOutSpec()
#    sim.genOutSpec('/Users/cdelacroix/Desktop/' + str(sim._outspec['seed']))

#import cProfile as profile
#profile.run('sim.SurveySimulation.run_sim()','output.qcg')

# %timeit(sim.run_sim(), sim.reset_sim())



# if True:
#     print 'IWA', OS.observingModes[0]['IWA']
#     print 'OWA', OS.observingModes[0]['OWA']
#     DRM = SS.DRM
#     for i in range(len(DRM)):
#         if any(DRM[i]['plan_inds']):
#             print i
#             print 'det_WA', DRM[i]['det_WA']
#             print 'det_status', DRM[i]['det_status']
#             print 'det_t', DRM[i]['det_time']
#             print 'det_SNR', DRM[i]['det_SNR']
#             print 'char_status', DRM[i]['char_status']
#             print 'char_time', DRM[i]['char_time']
#             print 'char_SNR', DRM[i]['char_SNR']


# pdb
# for j in range(100):
#     print j
#     sim.TimeKeeping.__init__(**sim.TimeKeeping._outspec)
#     sim.SimulatedUniverse.planTime = np.zeros(sim.SimulatedUniverse.nPlans)*u.d
#     sim.SurveySimulation.run_sim()


#def save_object(obj, filename):
#    with open(filename, 'wb') as output:
#        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
