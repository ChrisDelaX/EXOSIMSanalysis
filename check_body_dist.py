import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import copy
try:
    import cPickle as pickle
except:
    import pickle

scriptfile = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts/sampleScript_coron.json')
sim = EXOSIMS.MissionSim.MissionSim(scriptfile, forceStaticEphem=True)

TK = sim.TimeKeeping
OS = sim.OpticalSystem
TL = sim.TargetList
Obs = sim.Observatory
print 'forceStaticEphem = ', Obs.forceStaticEphem


# startTime = TK.missionStart
# startTime = Time(60634, format='mjd', scale='tai')
startTime = Time(64693, format='mjd', scale='tai')

t0 = startTime.mjd
t1 = t0 + 365.*2
dt = 2.  #365./8#
currentTimes = Time(np.arange(t0, t1+1, dt), format='mjd', scale='tai')

printvalues = False
bodies = ['Sun','Earth','Moon']#,'Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto']
d_sun = np.zeros(currentTimes.size)*u.AU
d_earth = np.zeros(currentTimes.size)*u.AU
d_sc = np.zeros(currentTimes.size)*u.AU
for body in bodies:
    for i,currentTime in enumerate(currentTimes):
        r_sun = Obs.solarSystem_body_position(currentTime, body).to('AU')
        r_earth = r_sun - Obs.solarSystem_body_position(currentTime, 'Earth').to('AU')
        r_sc = r_sun - Obs.orbit(currentTime).to('AU')
        d_sun[i] = np.linalg.norm(r_sun, axis=1)*u.AU
        d_earth[i] = np.linalg.norm(r_earth, axis=1)*u.AU
        d_sc[i] = np.linalg.norm(r_sc, axis=1)*u.AU
    if printvalues:
        print '%s-Sun distances (4 seasons): %s'%(body, np.round(d_sun, 3))
        print '%s-Earth distances (4 seasons): %s'%(body, np.round(d_earth, 3))
        print '%s-WFIRST distances (4 seasons): %s'%(body, np.round(d_sc, 3))
     
    if body == 'Sun':
        SunWFIRST = copy.copy(d_sc)
    if body == 'Earth':
        SunEarth = copy.copy(d_sun)
        EarthWFIRST = copy.copy(d_sc)
    if body == 'Moon':
        SunMoon = copy.copy(d_sun)
        EarthMoon = d_earth
        MoonWFIRST = copy.copy(d_sc)


# create figures
################

x = currentTimes.mjd - currentTimes.mjd[0]

fig, ax = plt.subplots(figsize = (12, 6))
plt.plot(x, SunEarth, label='Sun-Earth')
plt.plot(x, SunMoon, label='Sun-Moon')
plt.plot(x, SunWFIRST, label='Sun-WFIRST')
plt.xlim(x[0],x[-1])
plt.grid('on')
plt.xlabel('Time (days)')
plt.ylabel('Distance (AU)')
plt.title('Distance to Sun in AU')
plt.legend(fancybox=True, framealpha=.3, loc='best')
plt.savefig("Distance_to_Sun", dpi=300, transparent=True)

fig, ax = plt.subplots(figsize = (12, 6))
plt.plot(x, EarthMoon, label='Earth-Moon')
plt.plot(x, EarthWFIRST, label='Earth-WFIRST')
plt.xlim(x[0],x[-1])
plt.grid('on')
plt.xlabel('Time (days)')
plt.ylabel('Distance (AU)')
plt.title('Distance to Earth in AU')
plt.legend(fancybox=True, framealpha=.3, loc='best')
plt.savefig("Distance_to_Earth", dpi=300, transparent=True)

fig, ax = plt.subplots(figsize = (12, 6))
plt.plot(x, EarthWFIRST, label='WFIRST-Earth')
plt.plot(x, MoonWFIRST, label='WFIRST-Moon')
plt.xlim(x[0],x[-1])
plt.ylim(0,0.016)
plt.grid('on')
plt.xlabel('Time (days)')
plt.ylabel('Distance (AU)')
plt.title('Distance to WFIRST in AU')
plt.legend(fancybox=True, framealpha=.3, loc='best')
plt.savefig("Distance_to_WFIRST", dpi=300, transparent=True)

# Original data
###############

orbit_datapath = '/Users/cdelacroix/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Observatory/L2_halo_orbit_six_month.p'
halo = pickle.load(open(orbit_datapath, 'rb'))
t_halo = halo['t'][:,0]/(2*np.pi)*u.year
r_halo = halo['state'][:,0:3]*u.AU

fig, ax = plt.subplots(figsize = (12, 6))
plt.plot(t_halo*365., r_halo[:,0])
plt.ylim(1,1.012)
plt.grid('on')
plt.xlabel('Time (days)')
plt.ylabel('Distance (AU)')
plt.title('Original halo data')
plt.savefig("halo_data", dpi=300, transparent=True)


