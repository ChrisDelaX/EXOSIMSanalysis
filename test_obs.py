import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord


folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/EXOSIMS/gitexosims/EXOSIMS/Scripts'))
filename = 'mike_bottom.json'
scriptfile = os.path.join(folder,filename)
koangle = 45*u.deg
bodyname = 'Moon'
nb = 30
curt = 60634.

# TEST 1: Static ephemerides
############################
print '\nTEST 1: static ephemerides'
sim = EXOSIMS.MissionSim.MissionSim(scriptfile, forceStaticEphem=True)
Obs = sim.Observatory
TL = sim.TargetList

# 1 star, 1 time
sInds = 0
currentTime = Time(curt, format='mjd', scale='tai')
print '\n1 star, 1 time:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# n stars, 1 time
sInds = np.arange(nb)
print '\nn stars, 1 time:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# n stars, n times
currentTime = Time(np.array([curt]*nb), format='mjd', scale='tai')
print '\nn stars, n times:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# 1 star, n times
sInds = 0
print '\n1 star, n times:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape
print '\n'

# TEST 2: JPL ephemerides
############################
print '\nTEST 2: JPL ephemerides'
reload(EXOSIMS)
reload(EXOSIMS.MissionSim)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile, forceStaticEphem=False)
Obs = sim.Observatory
TL = sim.TargetList

# 1 star, 1 time
sInds = 0
currentTime = Time(curt, format='mjd', scale='tai')
print '\n1 star, 1 time:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# n stars, 1 time
sInds = np.arange(nb)
print '\nn stars, 1 time:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# n stars, n times
currentTime = Time(np.array([curt]*nb), format='mjd', scale='tai')
print '\nn stars, n times:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape

# 1 star, n times
sInds = 0
print '\n1 star, n times:'
print 'orbit', Obs.orbit(currentTime).shape
print 'keepout', Obs.keepout(TL, sInds, currentTime, koangle).shape
print 'starprop', Obs.starprop(TL, sInds, currentTime).shape
print 'moon-earth', Obs.moon_earth(currentTime).shape
print 'body pos', Obs.solarSystem_body_position(currentTime, bodyname).shape
