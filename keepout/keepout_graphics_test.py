#!/usr/bin/env python
#
# Make plots of keepout regions
#
# For usage, use the -h option.
#
# Examples:
# python keepout_graphics_test.py -s 64693 -l 2.0 -d 2 -m $HOME/Desktop/keepout_long.mp4 ./HabEx_4m_SS_20170411-keepout.json
# python keepout_graphics_test.py -s 60900 -l 0.3 -d 5 -m $HOME/Desktop/keepout.mp4 ./HabEx_4m_SS_20170411-keepout.json
# Where:
#  -s is the start-time of the movie (MJD)
#  -l is the length in years
#  -d is the delta-t between frames in days
#  -m is the name of the movie
#  [optionally, -f is the directory name for .png movie-frame output]
#  [optionally, -c is the directory name for cumulative keepout output]
#
# author:
#  Christian Delacroix, Oct 2016
# updates:
#  Michael Turmon, Oct 2016
#  Christian Delacroix, Mar 2017
#  Christian Delacroix, May 2017 (added occulter koAngleMax)
#  turmon may 2017: plot cumulative results
#  Delacroix June 8 update to coordinate frame selection

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
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D

# keeps xlabel from being chopped off
rcParams.update({'figure.autolayout': True})

# select coordinate frame: 
# 1 = heliocentric equatorial
# 2 = heliocentric ecliptic
coord_frame = 2

def make_graphics(args):

    ## set up arguments
    # movie output: False, or a path (file)
    out_movie = args.out_movie
    # frame output: False, or a path (*directory*, not file)
    out_frames = args.out_frames
    # cumulative output: False, or a path (*directory*, not file)
    out_cume = args.out_cume
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
    # make aspect ratio close to 2 wide x 1 tall for a lat/lon plot
    fig.set_size_inches(6,3.4)
    
    if out_movie:
        # movie hooked to "fig"
        #plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
        FFMpegWriter = animation.writers['ffmpeg']
        movie = FFMpegWriter(fps=15)
        movie.setup(fig, out_movie, 200)
    if out_frames:
        if not os.path.isdir(out_frames):
            os.makedirs(out_frames)
    if out_cume:
        if not os.path.isdir(out_cume):
            os.makedirs(out_cume)
    
    ########################################
    #
    # Set up loop
    # 
    ########################################
    
    endTime = startTime + nyears*365.0
    currentTimes = Time(np.arange(startTime, endTime+1, dt), format='mjd', scale='tai')
    
    # for progress indication
    system_t0 = time()
    system_dt = 10.0 # seconds
    
    ########################################
    #
    # Loop over times
    # 
    ########################################
    
    # get modules and initialize arrays
    Obs = sim.Observatory
    OS = sim.OpticalSystem
    TL = sim.TargetList
    mode = OS.observingModes[0]
    koMin = Obs.koAngleMin.to('deg').value
    koMax = Obs.koAngleMax.to('deg').value
    sInds = np.arange(TL.nStars)
    evergood = np.zeros(TL.nStars)
    
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
        kogood, r_body, r_targ, _ = Obs.keepout(TL, sInds, currentTime, mode, 
                returnExtra=True)
        evergood += kogood
        # strings for later titles
        currentTime.out_subfmt = 'date_hm' # suppress sec and ms on string output
        time_iso = currentTime.iso
        time_mjd = '%.1f' % currentTime.mjd
        
        # calculate equatorial coordinates of bright bodies and target stars
        targ = SkyCoord(SkyCoord(r_targ[:,0], r_targ[:,1], r_targ[:,2],
                representation='cartesian').represent_as('spherical'), 
                frame='heliocentrictrueecliptic', unit='deg,deg,pc')
        body = SkyCoord(SkyCoord(r_body[:,0,0], r_body[:,0,1], r_body[:,0,2],
                representation='cartesian').represent_as('spherical'), 
                frame='heliocentrictrueecliptic', unit='deg,deg,AU')
        
        # HELIOCENTRIC EQUATORIAL FRAME
        if coord_frame == 1:
            # get astropy SkyCoordinates
            targ = SkyCoord(SkyCoord(r_targ[:,0], r_targ[:,1], r_targ[:,2],
                    representation='cartesian').represent_as('spherical'), 
                    frame='icrs', unit='deg,deg,pc')
            body = SkyCoord(SkyCoord(r_body[:,0,0], r_body[:,0,1], r_body[:,0,2],
                    representation='cartesian').represent_as('spherical'), 
                    frame='icrs', unit='deg,deg,AU')
            # get RA and DEC
            x1 = targ.ra.to('deg').value[kogood]
            y1 = targ.dec.to('deg').value[kogood]
            x2 = targ.ra.to('deg').value[~kogood]
            y2 = targ.dec.to('deg').value[~kogood]
            x3 = body.ra.to('deg').value
            y3 = body.dec.to('deg').value
            movie_xlabel = 'equatorial RA'
            movie_ylabel = 'equatorial DEC'
        
        # HELIOCENTRIC ECLIPTIC FRAME
        elif coord_frame == 2:
            # transform from equatorial to ecliptic
            r_targ1 = Obs.equat2eclip(r_targ, currentTime)
            r_body1 = Obs.equat2eclip(r_body[:,0,:], currentTime)
            # get astropy SkyCoordinates
            targ = SkyCoord(SkyCoord(r_targ1[:,0], r_targ1[:,1], r_targ1[:,2],
                    representation='cartesian').represent_as('spherical'), 
                    frame='heliocentrictrueecliptic', unit='deg,deg,pc')
            body = SkyCoord(SkyCoord(r_body1[:,0], r_body1[:,1], r_body1[:,2],
                    representation='cartesian').represent_as('spherical'), 
                    frame='heliocentrictrueecliptic', unit='deg,deg,AU')
            # get LON and LAT
            x1 = targ.lon.to('deg').value[kogood]
            y1 = targ.lat.to('deg').value[kogood]
            x2 = targ.lon.to('deg').value[~kogood]
            y2 = targ.lat.to('deg').value[~kogood]
            x3 = body.lon.to('deg').value
            y3 = body.lat.to('deg').value
            movie_xlabel = 'ecliptic LON'
            movie_ylabel = 'ecliptic LAT'
        
        
        #import pdb;pdb.set_trace()
        # this is the earth-to-observatory distance
        print time_iso, body[2].distance.to('au')
        
        # planet names for labels
        #   now using L for Luna (moon) due to over-use of M
        names  = ['S',     'L', 'E', 'Mer', 'Ven', 'Mar', 'Jup', 'Sat', 'Ura', 'Nep', 'Plu']
        colors = ['orange','g', 'b', 'c',   'c',   'c',   'c',   'c',   'c',   'c',   'c']
        nBodies = np.size(names)
        koangles = np.ones(nBodies)*koMin
        koangles[3:] = 1. # smaller bodies angle is 1 deg
        
        if out_frames or out_movie:
            # clear previous figure
            plt.clf()
            # get axes
            ax = fig.gca()
            # faint line at equator/ecliptic
            plt.hlines(0, 0, 360, colors='gray', linestyles='dashed', linewidth=0.5)
            
            for j in range(nBodies):
                alpha = 0.2 if j <= 2 else 0.9
                if OS.haveOcculter and (j == 0):
                    width = 360-koMax
                    # change alpha to distinguish a bit
                    ring = patches.Wedge((x3[j],y3[j]), 360, 0, 360, width=width, 
                                                    color=colors[j], alpha=alpha/2)
                    ax.add_artist(ring)
                    ax.add_artist(patches.Circle((x3[j],y3[j]), koMax, 
                                                    color='gray', fill=False))
                    # overflow to right: replot on left
                    if x3[j] + koMax > 360:
                        ax.add_artist(patches.Circle((x3[j]-360,y3[j]), koMax, 
                                                    color='gray', fill=False))
                    # overflow to left: replot on right
                    if x3[j] - koMax < 0:
                        ax.add_artist(patches.Circle((x3[j]+360,y3[j]), koMax, 
                                                    color='gray', fill=False))
                disk = plt.Circle((x3[j],y3[j]), koangles[j], 
                                                    color=colors[j], alpha=alpha)
                ax.add_artist(disk)
                # disk underflows left edge of plot window: replot on right
                if x3[j] < koangles[j]:
                    disk = plt.Circle((x3[j] + 360.0, y3[j]), koangles[j], 
                                                    color=colors[j],alpha=alpha)
                    ax.add_artist(disk)
                # disk overflows right edge of plot window: replot on left
                if x3[j] > 360.0 - koangles[j]:
                    disk = plt.Circle((x3[j] - 360.0, y3[j]), koangles[j], 
                                                    color=colors[j],alpha=alpha)
                    ax.add_artist(disk)
                # marker for the text name
                ax.add_artist(plt.Text(text=names[j], x=x3[j], y=y3[j]))
                # add a marker to indicate the center for large objects
                if koangles[j] > 10:
                    ax.add_artist(plt.Text(text='+', x=x3[j], y=y3[j], 
                                                    color='gray',
                                                    horizontalalignment='center',
                                                    verticalalignment='center'))
                
                plt.minorticks_on()
                plt.xlim(0, 360)
                plt.ylim(-90, 90)
                plt.xlabel('%s [deg]' % movie_xlabel, fontsize=16)
                plt.ylabel('%s [deg]' % movie_ylabel, fontsize=16)
                day = currentTime.mjd - startTime
                plt.title('%s -- MJD %s -- Day #%.1f' % (time_iso, time_mjd, day))
                plt.tick_params(axis='both', which='major', labelsize=12)
                
                scatter_props = dict(s=10, edgecolors='none')
                plt.scatter(x1, y1, c=[0,.2,1], **scatter_props) # Visible = blue
                plt.scatter(x2, y2, c=[1,.2,0], **scatter_props) # Blocked = red
            
            # file output
            if out_frames:
                fn = os.path.join(out_frames, 'target_%04d.png' % i)
                plt.savefig(fn, dpi=300, transparent=True)
            # movie output
            if out_movie:
                movie.grab_frame()
    
    sys.stdout.write('\n\n')
    
    # need to finish the movie now, because it is linked to the current figure
    if out_movie:
        print 'Movie in %s' % movie.outfile
        movie.finish()
    if out_frames:
        print "Output image sequence in `%s'." % out_frames
    
    # cumulative outputs
    if out_cume:
        # all target coords
        xa = targ.icrs.ra.to('deg').value
        ya = targ.icrs.dec.to('deg').value
        # cumulative good observations
        cume_good = evergood * (1.0 / Ntime)
        cume_good_s = np.sort(cume_good)
        # times
        t_start = currentTimes[0]
        t_end   = currentTimes[-1]
        t_start.out_subfmt = 'date'
        t_end.out_subfmt = 'date'
        
        ## map-format plot of observability
        plt.clf()
        scatter_props = dict(s=20, cmap='plasma')
        plt.scatter(xa, ya, c=evergood/(1.0*Ntime), **scatter_props)
        # style it
        plt.minorticks_on()
        plt.xlim(0, 360);
        plt.ylim(-90, 90);
        plt.xlabel('RA [deg]', fontsize=16)
        plt.ylabel('DEC [deg]', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.title('Cumulative Obs. Map: %s -- %s' % (t_start.iso, t_end.iso))
        plt.colorbar()
        # save it
        fn = os.path.join(out_cume, 'cume_map.png')
        plt.savefig(fn, dpi=300, transparent=True)
        
        ## cumulative plot of observability
        plt.clf()
        plt.step(cume_good_s, np.linspace(0, 1, TL.nStars))
        plt.grid(True)
        # style it
        plt.xlabel('Portion of time observable', fontsize=16)
        plt.ylabel('Portion of targets observed', fontsize=16)
        plt.title('Cumulative Targets: %s -- %s' % (t_start.iso, t_end.iso))
        # 0..1 and 0..nStars y-scales
        ax1 = plt.gca()
        ax1.set_ylim((0.0, 1.0))
        #ax1.grid(True)
        par1 = ax1.twinx()
        par1.set_ylim((0, TL.nStars))
        # save it
        fn = os.path.join(out_cume, 'cume_hist.png')
        plt.savefig(fn, dpi=300, transparent=True)
        
        ## map-format plot of magnitudes
        # sort by mag, so that bright (=large) stars plot last/first
        idx = np.argsort(TL.Vmag)
        x = TL.coords.ra[idx]
        y = TL.coords.dec[idx]
        # sqrt() compresses a little
        z = np.sqrt(1e4*10**(-0.4*TL.Vmag[idx]))
        # make the plot
        plt.clf()
        # size is as above, color is Vmag directly
        #  this allows colorbar to work
        plt.scatter(x, y, c=TL.Vmag[idx], s=z, edgecolors='none')
        plt.colorbar()
        plt.minorticks_on()
        plt.xlim(0, 360);
        plt.ylim(-90, 90);
        # labels
        plt.xlabel('RA [deg]', fontsize=16)
        plt.ylabel('DEC [deg]', fontsize=16)
        plt.title('Target V Magnitude')
        plt.tick_params(axis='both', which='major', labelsize=14)
        fn = os.path.join(out_cume, 'magnitudes.png')
        plt.savefig(fn, dpi=300, transparent=True)
    
    #import pdb; pdb.set_trace();
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Make EXOSIMS keepout graphics.",
            epilog='At least one of -m or -f or -c should be given.')
    parser.add_argument('script', metavar='SCRIPT', help='json script')
    parser.add_argument('-c', '--cume', help='directory for output of cumulative info',
            dest='out_cume', metavar='DIR', default='')
    parser.add_argument('-f', '--frames', help='directory for output of still image frames',
            dest='out_frames', metavar='DIR', default='')
    parser.add_argument('-m', '--movie', help='file for output movie',
            dest='out_movie', metavar='FILE', default='')
    parser.add_argument('-s', '--start', help='start time [mjd], default = %(default).0f [20-Nov-2024]',
            type=float, dest='t0', default=60634.0)
    parser.add_argument('-l', '--length', help='length [years], default = %(default).2f',
            type=float, dest='years', default=0.5)
    parser.add_argument('-d', '--delta', help='delta-time [days], default = %(default).2f',
            type=float, dest='delta_t', default=5.0)
    
    args = parser.parse_args()
    
    if not args.out_movie and not args.out_frames and not args.out_cume:
        print 'Reminder: Need one of -m or -f or -c in order to produce output.'
        # sys.exit(1)
    
    make_graphics(args)
