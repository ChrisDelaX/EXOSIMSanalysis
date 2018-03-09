import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
import astropy.constants as const
from EXOSIMS.util.keplerSTM import planSys

filename = 'sampleScript_coron.json'
scriptfolder = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts')
scriptfile = os.path.join(scriptfolder, filename)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)

SU = sim.SimulatedUniverse
TL = sim.TargetList

sInd = 0
dt = 1e5*u.d
nbmax = 1e2
for attempt in range(int(nbmax)):
    pInds = np.where(SU.plan2star == sInd)[0]
    nPlans = pInds.size
    print 'attempt number %s with a %s-planet system (%s-planet universe)'%(attempt, 
            nPlans, SU.nPlans)
    try:
        SU.propag_system(sInd, dt)
    except:
        Ms = TL.MsTrue[[sInd]]
        Mp = SU.Mp[pInds]
        mu = (const.G*(Mp + Ms)).to('AU3/day2').value
        r0 = SU.r[pInds].to('AU').value
        v0 = SU.v[pInds].to('AU/day').value
        x0 = np.reshape(np.concatenate((r0, v0), axis=1), nPlans*6)
        print 'plansys error occured: mu and x0 have been saved'
        break
    else:
        # if plansys error doesn't occur, reset_sim and try again
        sim.reset_sim()
else:
    raise ValueError("plansys error didn't occur after %s attempts"%nbmax)

