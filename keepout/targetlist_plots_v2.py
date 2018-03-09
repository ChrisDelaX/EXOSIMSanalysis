import sys
import os
import os.path
from time import time
import EXOSIMS
import EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.animation as animation


########################################
#
# Filenames
#
########################################

# Script

#folder = os.path.normpath(os.path.expandvars('$HOME/CODE/Python/exomissionsim/EXOSIMS/Scripts'))
script_folder = '.'
script_filename = 'sampleScript_coron_v2.json'
script_path = os.path.join(script_folder, script_filename)

# load script
sim = EXOSIMS.MissionSim.MissionSim(script_path)

# Graphics -- need to associate the movie with the figure now
plt.close('all')
fig = plt.gcf()

# movie output: False, or a path (file)
#out_movie = False
out_movie = '/tmp/keepout_test.mp4'
# file output: False, or a path (*directory*, not file)
out_files = False
#out_files = '/tmp/keepout'

if out_movie:
    # movie hooked to "fig"
    #plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
    FFMpegWriter = animation.writers['ffmpeg']
    movie = FFMpegWriter(fps=15)
    movie.setup(fig, out_movie, 200)
if out_files:
    if not os.path.isdir(out_files):
        os.makedirs(out_files)


########################################
#
# Set up loop
# 
########################################

nyears = 0.5
startTime = 60634.0 + 200.0 # mjd
dt = 1.0 # days
endTime = startTime + nyears*365.0
currentTimes = Time(np.arange(startTime, endTime, dt), format='mjd', scale='tai')

# get objects
OS = sim.OpticalSystem
Obs = sim.Observatory
TL = sim.TargetList
sInds = np.arange(TL.nStars)
koangle = OS.telescopeKeepout
evergood = np.zeros(TL.nStars)

# for progress indication
system_t0 = time()
system_dt = 10.0 # seconds

########################################
#
# Loop over times
# 
########################################

print 'Beginning keepout computation.'
Ntime = len(currentTimes)
for i in range(Ntime):
    # print progress if needed
    if time() > system_t0 + system_dt:
        sys.stdout.write('\n[%3d/%3d]' % (i+1, Ntime))
        system_t0 = time()
    else:
        sys.stdout.write('.')
        sys.stdout.flush()

    currentTime = currentTimes[i]
    currentTime = currentTime.reshape(currentTime.size)

    # strings for later titles
    currentTime.out_subfmt = 'date_hm' # suppress sec and ms on string output
    time_iso = currentTime[0].iso
    time_mjd = '%.1f' % currentTime[0].mjd

    r_sc = Obs.orbit(currentTime)
    kogood = Obs.keepout(TL, sInds, currentTime, r_sc, koangle)
    evergood += kogood
    
    r_targ = Obs.starprop(TL, sInds, currentTime)
    NewCoords = SkyCoord(r_targ[:,0],r_targ[:,1],r_targ[:,2],representation='cartesian')
    NewCoords = NewCoords.heliocentrictrueecliptic.icrs
    # visible targets, and kept-out targets
    x1 = NewCoords.ra [ kogood]
    y1 = NewCoords.dec[ kogood]
    x2 = NewCoords.ra [~kogood]
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
    
    # planet names for labels
    names = ['S', 'M', 'E',
        'Mer', 'Ven', 'Mar',
        'Jup', 'Sat', 'Ura', 'Nep', 'Plu']

    nBodies = r_bright.shape[0]
    koangles = np.ones(nBodies)*koangle
    koangles[3:] = 1*u.deg
    colors = ['orange','g','b','c','c','c','c','c','c','c','c']
    
    body = SkyCoord(r_bright[:,0,0],
                    r_bright[:,0,1],
                    r_bright[:,0,2],
                    representation='cartesian')
    body = body.heliocentrictrueecliptic.icrs
    
    #import pdb; pdb.set_trace()

    plt.clf()
    ax = fig.gca()
    for j in range(nBodies):
        alpha = 0.2 if j <= 2 else 0.9
        disk = plt.Circle((body.ra[j].value, body.dec[j].value),
                    koangles[j].value, color=colors[j],alpha=alpha)
        ax.add_artist(disk)
        ax.add_artist(plt.Text(text=names[j], x=body.ra[j].value, y=body.dec[j].value))
        # disk underflows left edge of plot window: replot on right
        if body.ra[j].value < koangles[j].value:
            disk = plt.Circle((body.ra[j].value + 360.0, body.dec[j].value), 
                    koangles[j].value, color=colors[j],alpha=alpha)
            ax.add_artist(disk)
        # disk overflows right edge of plot window: replot on left
        if body.ra[j].value > 360.0 - koangles[j].value:
            disk = plt.Circle((body.ra[j].value - 360.0, body.dec[j].value), 
                    koangles[j].value, color=colors[j],alpha=alpha)
            ax.add_artist(disk)
    
    plt.minorticks_on()
    plt.xlim(0, 360);
    plt.ylim(-90, 90);
    plt.xlabel('RA [deg]', fontsize=16)
    plt.ylabel('DEC [deg]', fontsize=16)
    day = currentTime.value - startTime
    plt.title('%s -- MJD %s -- Day #%.1f' % (time_iso, time_mjd, day[0]))
    plt.tick_params(axis='both', which='major', labelsize=14)
    
    scatter_props = dict(s=10, edgecolors='none')
    plt.scatter(x1, y1, c=[0,.2,1], **scatter_props) # Visible = blue
    plt.scatter(x2, y2, c=[1,.2,0], **scatter_props) # Blocked = red
    # file output
    if out_files:
        fn = os.path.join(out_files, 'target_%04d.png' % i)
        plt.savefig(fn, dpi=300, transparent=True)
    # movie output
    if out_movie:
        movie.grab_frame()

sys.stdout.write('\n\n')

if out_movie:
    print 'Movie in %s' % movie.outfile
    movie.finish()

if out_files:
    print 'Output image sequence in %s' % out_files

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

