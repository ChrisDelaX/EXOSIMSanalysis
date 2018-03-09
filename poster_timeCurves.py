import os.path,EXOSIMS,EXOSIMS.MissionSim
#scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','HabX_4m_HZ_20160621_1114.json')
#scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','template_WFIRST_EarthTwinHabZone.json')
#scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','template_WFIRST_KnownRV.json')
#scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','sampleScript_coron.json')
folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/exomissionsim/MyFiles/Scripts/'))
scriptfile = 
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)

#import cProfile as profile
#profile.run('sim.SurveySimulation.run_sim()','output.qcg')
%pylab

import os, inspect, copy, numbers, json, pandas, h5py 
import numpy as np
import cPickle as pickle
import astropy.io.fits as fits
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.io.votable import parse
from scipy.interpolate import interp1d, interp2d, griddata
from EXOSIMS.util.deltaMag import deltaMag
from EXOSIMS.util.get_module import get_module
from EXOSIMS.util.eccanom import eccanom
from astropy.coordinates import SkyCoord
Obs = sim.Observatory
OS = sim.OpticalSystem
TL = sim.TargetList
TK = sim.TimeKeeping
SS = sim.SurveySimulation
SU = sim.SimulatedUniverse
PP = sim.PostProcessing
PPop = sim.PlanetPopulation
PPMod = sim.PlanetPhysicalModel
BS = sim.BackgroundSources
ZL = sim.ZodiacalLight
specs = json.loads(open(scriptfile).read())
SC = get_module(specs['modules']['StarCatalog'],'StarCatalog')(**specs)


# 
# # create set of planets with random planet delta magnitudes
# IWA = OS.IWA.value
# OWA = OS.OWA.value
# pInds = np.array(range(SU.nPlans))
# allWA = SU.get_current_WA(pInds)
# WA = allWA[allWA>OS.IWA]
# I = SU.I[pInds][allWA>OS.IWA]
# sInd = SU.plan2star[allWA>OS.IWA]
# I = np.array([y for (x,y) in sorted(zip(WA,I.value))])*u.deg
# sInd = np.array([y for (x,y) in sorted(zip(WA,sInd))])
# WA.sort()
# dMag = np.random.normal(21.3,1,WA.size); # Bijan's values
# dMag.sort()
# 
# # different cameras
# cam = np.array(['CCD','AM-EMCCD','PC-EMCCD','sCMOS']);
# ENF = np.array([1., np.sqrt(2), 1., 2.]);
# q_cic = np.array([1e-3, 1e-3, 1e-3, 0.]);
# read_noise = np.array([3., 16./500, 16./500, 1.]);
# 
# %pylab
# lwz = 3;
# fsz = 12;
# grid('on');xlabel('planet-star dMag',fontsize=1.5*fsz);ylabel('charact time (days) - SNR=5',fontsize=1.5*fsz)
# xmax = 25 #t.size-1
# #plot([1,xmax],[30,30],'k--',linewidth=lwz)
# yscale('log');ylim(6e-2,100);
# xlim(20,22)
# #xticks([0,100,200,365], fontsize = fsz)
# #yticks([1e-1,1e0,1e1], fontsize = fsz)
# PP.SNchar = 5
# for i in range(3):
#     OS.scienceInstruments[1]['ENF'] = ENF[i];
#     OS.scienceInstruments[1]['CIC'] = q_cic[i];
#     OS.scienceInstruments[1]['sread'] = read_noise[i];
#     tt = OS.calc_charTime(TL,sInd,dMag,0,0/u.arcsec**2,0/u.arcsec**2).value;
#     t = tt[tt>0]; 
#     t.sort();
#     x = range(len(t)+1)[1:];
#     y = np.empty(len(t))
#     for j in range(len(t)):
#         if j == 0:
#             y[j]=t[j];
#         else:
#             y[j]=t[j]+y[j-1];
#     plot(dMag,t,'-',label=cam[i],linewidth=lwz, markersize=3*lwz);
# #print t
# legend(fancybox=True, framealpha=0.3, loc='lower right', fontsize=fsz);
# tight_layout()
# savefig('MyFiles/'+str(int(OS.Spectro['lam'].value))+'.png', dpi=300, transparent=True);
# 

