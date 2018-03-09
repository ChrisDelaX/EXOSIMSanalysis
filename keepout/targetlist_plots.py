import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
%pylab

folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/exomissionsim/EXOSIMS/Scripts'))
filename = 'sampleScript_coron.json'
scriptfile = os.path.join(folder,filename)
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)

OS = sim.OpticalSystem
Obs = sim.Observatory
TL = sim.TargetList

nyears = 3.0
startTime = 60634.0
dt = 5.0
endTime = startTime + nyears*365.0
currentTimes = Time(np.arange(startTime, endTime, dt), format='mjd', scale='tai')

sInds = np.arange(TL.nStars)
koangle = OS.telescopeKeepout

evergood = np.zeros(TL.nStars)
for i in range(len(currentTimes)):
    currentTime = currentTimes[i]
    currentTime = currentTime.reshape(currentTime.size)
    r_sc = Obs.orbit(currentTime)
    kogood = Obs.keepout(TL, sInds, currentTime, r_sc, koangle)
    evergood += kogood
    
    # find target coordinates in the cartesian frame
    r_targ = Obs.starprop(TL, sInds, currentTime)
    NewCoords = SkyCoord(r_targ[:,0],r_targ[:,1],r_targ[:,2],representation='cartesian')
    # transform to ecliptic frame
    NewCoords = NewCoords.heliocentrictrueecliptic
    # transform to equatorial frame
    NewCoords = NewCoords.icrs
    
    x1 = NewCoords.ra[kogood]
    y1 = NewCoords.dec[kogood]
    x2 = NewCoords.ra[~kogood]
    y2 = NewCoords.dec[~kogood]
    
    r_bright = np.array([np.zeros(r_sc.shape), # sun
        Obs.solarSystem_body_position(currentTime, 'Moon').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Earth').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Mercury').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Venus').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Mars').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Jupiter').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Saturn').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Uranus').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Neptune').T.to('km').value,
        Obs.solarSystem_body_position(currentTime, 'Pluto').T.to('km').value])*u.km
    r_bright -= r_sc
    
    nBodies = r_bright.shape[0]
    koangles = np.ones(nBodies)*koangle
    koangles[3:] = 1*u.deg
    colors = ['orange','g','b','c','c','c','c','c','c','c','c']
    
    # find body coordinates in the cartesian frame
    body = SkyCoord(r_bright[:,0,0],r_bright[:,0,1],r_bright[:,0,2],representation='cartesian')
    # transform to ecliptic frame
    body = body.heliocentrictrueecliptic
    # transform to equatorial frame
    body = body.icrs
    
    close('all')
    fig = gcf()
    ax = fig.gca()
    for j in range(nBodies):
        alpha = 0.2 if j <= 2 else 0.9
        disk = Circle((body.ra[j].value, body.dec[j].value), \
                koangles[j].value, color=colors[j],alpha=alpha)
        ax.add_artist(disk)
        if body.ra[j].value < koangles[j].value:
            disk = Circle((body.ra[j].value + 360.0, body.dec[j].value), \
                    koangles[j].value, color=colors[j],alpha=0.2)
            ax.add_artist(disk)
        if body.ra[j].value > 360.0 - koangles[j].value:
            disk = Circle((body.ra[j].value - 360.0, body.dec[j].value), \
                    koangles[j].value, color=colors[j],alpha=0.2)
            ax.add_artist(disk)
    
    minorticks_on()
    xlim(0,360);ylim(-90,90);
    xlabel('RA (deg)', fontsize=16)
    ylabel('DEC (deg)', fontsize=16)
    day = currentTime.value - startTime
    title('day = '+str(day))
    tick_params(axis='both', which='major', labelsize=14)
    
    sz = 10
    scatter(x1,y1,c=[0,.2,1],s=sz,edgecolors='none') # green
    scatter(x2,y2,c=[1,.2,0],s=sz,edgecolors='none') # red
    savefig('target'+str(i)+'.png', dpi=300, transparent=True);




##############

# # Plot of magnitudes:
#
# idx = np.argsort(TL.Vmag)
# x = TL.coords.ra[idx]
# y = TL.coords.dec[idx]
# z = 1e4*10**(-0.4*TL.Vmag[idx])
# 
# scatter(x,y,c=z,s=z, edgecolors='none')
# 
# minorticks_on()
# xlim(0,360);ylim(-90,90);
# xlabel('RA (deg)', fontsize=16)
# ylabel('DEC (deg)', fontsize=16)
# tick_params(axis='both', which='major', labelsize=14)
# 
# savefig("target1.png", dpi=300, transparent=True);

