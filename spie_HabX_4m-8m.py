%pylab

import os, inspect, copy, numbers, json, pandas 
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
from scipy.stats import norm

path = "MyFiles/HabX_4m/"

OBSs = []
DETs = []
CHARs = []
UNs = []
cnt = 0
for name in (os.listdir(path)):
    if name not in '.DS_Store':
        cnt += 1
        scriptfile = path + name
        specs = json.loads(open(scriptfile).read())
        DRM = specs['DRM']
        nobs = len(DRM)
        OBSs.append(nobs)
        ds = np.zeros(nobs)
        for i in range(nobs):
            if DRM[i]['det_status'] > 0:
                ds[i] = DRM[i]['det_status'][0]
        ds = np.array(ds)
        det = np.array([int(ds[i] >0) for i in range(nobs)])
        ndet = np.sum(det)
        DETs.append(ndet)
        cs = [DRM[i].get('char_1_time',np.nan) for i in range(nobs)]
        char = np.array([int(cs[i] >0) for i in range(nobs)])
        nchar = np.sum(char)
        CHARs.append(nchar)
        pInds = np.zeros(nobs)
        for i in range(nobs):
            if DRM[i]['plan_inds'] != []:
                pInds[i] = DRM[i]['plan_inds'][0]
        pInds = np.array(pInds)
        UNs.append(size(pInds[ds!=0]))

print np.mean(OBSs), np.mean(DETs), np.mean(CHARs), np.mean(UNs)


path = "MyFiles/HabX_8m/"

OBS2s = []
DET2s = []
CHAR2s = []
UN2s = []
for name in (os.listdir(path)):
    if name not in '.DS_Store':
        scriptfile = path + name
        specs = json.loads(open(scriptfile).read())
        DRM = specs['DRM']
        nobs = len(DRM)
        OBS2s.append(nobs)
        ds = [DRM[i]['det_status'] for i in range(nobs)]
        det = np.array([int(ds[i] >0) for i in range(nobs)])
        ndet = np.sum(det)
        DET2s.append(ndet)
        cs = [DRM[i].get('char_1_time',np.nan) for i in range(nobs)]
        char = np.array([int(cs[i] >0) for i in range(nobs)])
        nchar = np.sum(char)
        CHAR2s.append(nchar)
        pInds = [DRM[i]['plan_inds'] for i in range(len(DRM))]
        vp = []
        wp = []
        for v in pInds:
            if v:
                vp.append(v) 
            for w in v:
               wp.append(w) 
        UN2s.append(size(np.unique(wp)))


print np.mean(OBS2s), np.mean(DET2s), np.mean(CHAR2s), np.mean(UN2s)

#D1 = UNs
#D2 = UN2s
#xmax = 11
D1 = DETs
D2 = CHARs
xmax = 18
bins = range(xmax+1)
(mu, sigma) = norm.fit(D1)
y1 = mlab.normpdf( bins, mu, sigma)
(mu, sigma) = norm.fit(D2)
y2 = mlab.normpdf( bins, mu, sigma)

close('all')
figure(2)
grid('on')
fsz = 18
plot(bins, y1, 'b-o', linewidth=2, markersize=5, label='KasdinBraems')
plot(bins, y2, 'r-o', linewidth=2, markersize=5, label='Nemati')
#xlabel('Unique planet detections', fontsize=fsz)
xlabel('Total planet detections', fontsize=fsz)
ylabel('Normalized frequency', fontsize=fsz)
xlim(0,xmax)
tick_params(axis='both', which='major', labelsize=fsz)
legend(fancybox=True, framealpha=0.3, loc='upper right', fontsize=fsz);

#savefig("unique.png", dpi=300, transparent=True);
savefig("detections.png", dpi=300, transparent=True);

