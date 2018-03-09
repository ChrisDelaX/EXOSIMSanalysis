import sys,os.path
import numpy as np
#sys.path.insert(0,"/home/dean/Documents/exosims/EXOSIMS-master/EXOSIMS")
sys.path.append("/home/dean/Documents/exosims/EXOSIMS/EXOSIMS")
sys.path
import EXOSIMS,EXOSIMS.MissionSim
import astropy.units as u
from EXOSIMS.util.get_module import get_module
#folder = os.path.normpath(os.path.expandvars('$HOME/Documents/exosims/EXOSIMS/EXOSIMS/Scripts'))
folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/exomissionsim/MyFiles/Scripts'))
#filename = 'sampleScript_coron.json'
#filename= 'template_WFIRST_KnownRV_mod.json'
filename = 'sampleScript_integrationTime.json'
scriptfile = os.path.join(folder,filename)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
#sim.SurveySimulation.run_sim()


# Choose observing modes selected for detection (default marked with a flag),
#detMode = filter(lambda mode: mode['detectionMode'] == True, OS.observingModes)[0]
#sInd = None


#possible target list
#SU = get_module(specs['modules'] \ ['SimulatedUniverse'],'SimulatedUniverse')(**specs)

#SU = get_module(specs['modules']\
#                ['SimulatedUniverse'],'SimulatedUniverse')(**specs)
TL = sim.TargetList
#TL = self.TargetList
#t_dets = TL.nStars
#print t_dets


sInds = np.arange(TL.nStars)
print "sInds is "
print sInds


OS = sim.OpticalSystem
ZL = sim.ZodiacalLight
Obs = sim.Observatory
TK = sim.TimeKeeping
startTime = TK.currentTimeAbs+np.zeros(len(sInds))*u.d
r_sc = Obs.orbit(startTime)


mode = filter(lambda mode: mode['detectionMode'] == True, OS.observingModes)[0]
fZ = ZL.fZ(TL, sInds, mode['lam'], r_sc[sInds])
fEZ = ZL.fEZ0
t_dets= OS.calc_maxintTime(TL,sInds,fZ,fEZ,mode)



#Redo calc but with startTime as a single value
startTime = TK.currentTimeAbs


r_sc = Obs.orbit(startTime)
fZ = ZL.fZ(TL, sInds, mode['lam'], r_sc)
t_dets= OS.calc_maxintTime(TL,sInds,fZ,fEZ,mode)
print t_dets

#for i = 1:length target list 
# Acquire the NEXT TARGET star index and create DRM
#DRM, sInd, t_det = sim.next_target(sInd, detMode)
#t_det[i] = 
#end for


#spectroModes = filter(lambda mode: 'spec' in mode['inst']['name'], OS.observingModes)
#	if np.any(spectroModes):
#            charMode = spectroModes[0]
#        # if no spectro mode, default char mode is first observing mode
#        else:
#            charMode = OS.observingModes[0]
#characterized, charSNR, t_char = self.observation_characterization(sInd, charMode)

#sim.OpticalSystem.calc_intTime()
