import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

scriptfile = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/Scripts/HabEx_4m_SAG13_20170610.json')
plt.close('all')

sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
PPop = sim.PlanetPopulation
SU = sim.SimulatedUniverse

coeffs = PPop.SAG13coeffs.reshape(PPop.SAG13coeffs.size,order='F').tolist()

###########################################################
# FIRST generate a universe with default sma range values, 
# corresponding to period=[10,640] day
###########################################################
fig, ax = plt.subplots(figsize=(6,6))
myAxes = [PPop.Trange.to('day').value[0], PPop.Trange.to('day').value[1],
        PPop.Rprange.to('earthRad').value[0], PPop.Rprange.to('earthRad').value[1]]
im = ax.imshow(PPop.eta,origin='lower', extent=myAxes, aspect='auto')
ax.loglog()
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
plt.xticks([10,100,1000])
plt.yticks([1,2,5,10])
plt.xlabel('period T (day)')
plt.ylabel('radius (earth rad)')
fig.colorbar(im)
plt.title('coeffs: %s'%coeffs)
plt.savefig('SAG13_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)

fig, ax = plt.subplots(figsize=(6,6))
plt.hist(SU.Rp, bins=np.logspace(min(np.log(SU.Rp.value)), max(np.log(SU.Rp.value)), 50))
plt.xscale('log')
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
plt.xlim(0.5,20)
plt.xticks([0.5,1,2,5,10,20])
plt.xlabel('radius Rp (earth rad)')
plt.ylabel('# planets')
plt.savefig('radius_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)

fig, ax = plt.subplots(figsize=(6,6))
plt.hist(SU.T, bins=np.logspace(min(np.log(SU.T.value)), max(np.log(SU.T.value)), 50))
plt.xscale('log')
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
plt.xlim(0.02,4)
plt.xticks([0.02,0.1,0.3,1,4])
plt.xlabel('period T (year)')
plt.ylabel('# planets')
plt.savefig('period_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)

###########################################################
# THEN, run a new simulation with specific sma range values
###########################################################
reload(EXOSIMS)
reload(EXOSIMS.MissionSim)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile, arange=[0.2, 50])
PPop = sim.PlanetPopulation
SU = sim.SimulatedUniverse

fig, ax = plt.subplots(figsize=(6,6))
myAxes = [PPop.Trange.to('day').value[0], PPop.Trange.to('day').value[1],
        PPop.Rprange.to('earthRad').value[0], PPop.Rprange.to('earthRad').value[1]]
im = ax.imshow(PPop.eta, origin='lower', extent=myAxes, aspect='auto')
ax.loglog()
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
plt.xticks([350,380,410,440,470])
plt.yticks([1,2,5,10])
plt.xlabel('period T (day)')
plt.ylabel('radius (earth rad)')
fig.colorbar(im)
plt.title('coeffs: %s'%coeffs)
plt.savefig('SAG13_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)

fig, ax = plt.subplots(figsize=(6,6))
plt.hist(SU.Rp, bins=np.logspace(min(np.log(SU.Rp.value)), max(np.log(SU.Rp.value)), 50))
plt.xscale('log')
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
plt.xlim(0.5,20)
plt.xticks([0.5,1,2,5,10,20])
plt.xlabel('radius Rp (earth rad)')
plt.ylabel('# planets')
plt.savefig('radius_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)

fig, ax = plt.subplots(figsize=(6,6))
plt.hist(SU.T, bins=np.logspace(min(np.log(SU.T.value)), max(np.log(SU.T.value)), 50))
plt.xscale('log')
for axis in [ax.xaxis, ax.yaxis]: axis.set_major_formatter(ScalarFormatter())
#plt.xlim(0.9,1.3)
#plt.xticks([0.9,1,1.1,1.2,1.3])
plt.xlabel('period T (year)')
plt.ylabel('# planets')
plt.savefig('period_sma=%s.png'%PPop.arange.value.round(2), dpi=300, transparent=True)
