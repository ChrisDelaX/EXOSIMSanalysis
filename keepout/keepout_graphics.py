#!/usr/bin/env python
#
# Make plots of keepout regions
#
# For usage, use the -h option.
#
# Typical:
# python keepout_graphics.py -s 60900 -l 0.2 -d 0.5 -m $HOME/Desktop/keepout.mp4 ./sampleScript_coron_v2.json 
# Where:
# -s is the start-time of the movie (MJD)
# -l is the length in years
# -d is the delta-t between frames in days
# -m is the name of the movie
# [optionally, -f is the directory name for .png output]
#
# author:
#  Christian Delacroix, Oct 2016
# updates:
#  Michael Turmon, Oct 2016
#  Christian Delacroix, Mar 2017
#  Christian Delacroix, May 2017 (added occulter koAngleMax)
# 


import sys
import argparse
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
import matplotlib.patches as patches

def make_graphics(args):

    ## set up arguments
    # movie output: False, or a path (file)
    out_movie = args.out_movie
    # file output: False, or a path (*directory*, not file)
    out_files = args.out_files
    # times
    nyears = args.years
    dt = args.delta_t
    startTime = args.t0
    # script
    script_path = args.script
    
    ########################################
    #
    # Filenames
    #
    ########################################
    
    # load script
    sim = EXOSIMS.MissionSim.MissionSim(script_path)
    
    # Graphics -- need to associate the movie with the figure now
    plt.close('all')
    fig = plt.gcf()
    
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
    
    endTime = startTime + nyears*365.0
    currentTimes = Time(np.arange(startTime, endTime, dt), format='mjd', scale='tai')
    
    # get objects
    OS = sim.OpticalSystem
    Obs = sim.Observatory
    TL = sim.TargetList
    sInds = np.arange(TL.nStars)
    evergood = np.zeros(TL.nStars)
    koangle = Obs.koAngleMin
    
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
        
        # Build 'keepout-good' array
        currentTime = currentTimes[i]
        kogood, r_body, r_targ, _ = Obs.keepout(TL, sInds, currentTime, koangle, returnExtra=True)
        evergood += kogood
        # strings for later titles
        currentTime.out_subfmt = 'date_hm' # suppress sec and ms on string output
        time_iso = currentTime.iso
        time_mjd = '%.1f' % currentTime.mjd
        
        # Calculate (ra,dec) coordinates of visible targets, kept-out targets, and bright bodies
        targ = SkyCoord(r_targ[:,0],r_targ[:,1],r_targ[:,2],representation='cartesian')
        targ = targ.heliocentrictrueecliptic.icrs
        x1 = targ.ra.to('deg').value[kogood]
        y1 = targ.dec.to('deg').value[kogood]
        x2 = targ.ra.to('deg').value[~kogood]
        y2 = targ.dec.to('deg').value[~kogood]
        body = SkyCoord(r_body[:,0,0],r_body[:,0,1],r_body[:,0,2],representation='cartesian')
        body = body.heliocentrictrueecliptic.icrs
        x3 = body.ra.to('deg').value
        y3 = body.dec.to('deg').value
        
        # planet names for labels
        names = ['S','M','E','Mer','Ven','Mar','Jup', 'Sat', 'Ura', 'Nep', 'Plu']
        colors = ['orange','g','b','c','c','c','c','c','c','c','c']
        nBodies = np.size(names)
        koangles = np.ones(nBodies)*koangle.to('deg').value
        koangles[3:] = 1. # smaller bodies angle is 1 deg
        
        # create figure
        plt.clf()
        ax = fig.gca()
        for j in range(nBodies):
            alpha = 0.2 if j <= 2 else 0.9
            if OS.haveOcculter and (j==0):
                width = 360-Obs.koAngleMax.to('deg').value
                ring = patches.Wedge((x3[j],y3[j]), 360, 0, 360, width=width, color=colors[j], alpha=alpha)
                ax.add_artist(ring)
            disk = plt.Circle((x3[j],y3[j]), koangles[j], color=colors[j], alpha=alpha)
            ax.add_artist(disk)
            ax.add_artist(plt.Text(text=names[j], x=x3[j], y=y3[j]))
            # disk underflows left edge of plot window: replot on right
            if x3[j] < koangles[j]:
                disk = plt.Circle((x3[j] + 360.0, y3[j]), koangles[j], color=colors[j],alpha=alpha)
                ax.add_artist(disk)
            # disk overflows right edge of plot window: replot on left
            if x3[j] > 360.0 - koangles[j]:
                disk = plt.Circle((x3[j] - 360.0, y3[j]), koangles[j], color=colors[j],alpha=alpha)
                ax.add_artist(disk)
        
        plt.minorticks_on()
        plt.xlim(0, 360);
        plt.ylim(-90, 90);
        plt.xlabel('RA [deg]', fontsize=16)
        plt.ylabel('DEC [deg]', fontsize=16)
        day = currentTime.mjd - startTime
        plt.title('%s -- MJD %s -- Day #%.1f' % (time_iso, time_mjd, day))
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Make EXOSIMS keepout graphics.",
                                     epilog='At least one of -m or -f must be given.')
    parser.add_argument('script', metavar='SCRIPT', help='json script')
    parser.add_argument('-f', '--files', help='directory for output still images',
                      dest='out_files', metavar='DIR', default='')
    parser.add_argument('-m', '--movie', help='file for output movie',
                      dest='out_movie', metavar='FILE', default='')
    parser.add_argument('-s', '--start', help='start time [mjd], default = %(default).0f [20-Nov-2024]',
                      type=float, dest='t0', default=60634.0)
    parser.add_argument('-l', '--length', help='length [years], default = %(default).2f',
                      type=float, dest='years', default=0.5)
    parser.add_argument('-d', '--delta', help='delta-time [days], default = %(default).2f',
                      type=float, dest='delta_t', default=5.0)
    
    args = parser.parse_args()
    
    if not args.out_movie and not args.out_files:
        print 'Need one of -m or -f in order to produce output.'
        sys.exit(1)
    
    make_graphics(args)

##############

#     # Plot of magnitudes:
#     idx = np.argsort(TL.Vmag)
#     x = TL.coords.ra[idx]
#     y = TL.coords.dec[idx]
#     z = 1e4*10**(-0.4*TL.Vmag[idx])
#     plt.scatter(x,y,c=z,s=z, edgecolors='none')
#     plt.minorticks_on()
#     plt.xlim(0,360);plt.ylim(-90,90);
#     plt.xlabel('RA (deg)', fontsize=16)
#     plt.ylabel('DEC (deg)', fontsize=16)
#     plt.tick_params(axis='both', which='major', labelsize=14)
#     plt.savefig(os.path.normpath(os.path.expandvars("$HOME/Desktop/magnitudes.png")), dpi=300, transparent=True);

