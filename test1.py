# solar longitudes and latitudes of targets
lon = TL.coords.barycentrictrueecliptic.lon[sInds].value
lat = TL.coords.barycentrictrueecliptic.lat[sInds].value

#####

OS.calc_intTime(TL, range(TL.nStars), 22.5, (OS.IWA+OS.OWA)/2., 0./u.arcsec**2, 0./u.arcsec**2)

#####
c = TL.coords[0]

dec_rad = np.radians(c.dec).value
ra_rad = np.radians(c.ra).value
eps_rad = np.radians(23.439)
sinbeta = np.sin(dec_rad) * np.cos(eps_rad) - np.cos(dec_rad) * np.sin(eps_rad) * np.sin(ra_rad)
np.degrees(np.arcsin(sinbeta))

c.heliocentrictrueecliptic.lat.value

#########################

# print a dictionary:
def testfun(namedarg =  None, **specs):
    print 'namedarg : %s' % namedarg
    for key in specs.keys():
        print '%s : %s'%(key, specs[key])

specs = {}
specs['namedarg'] = 1
namedarg = 2
testfun(**specs)
testfun(namedarg=namedarg, **specs)


#########################

# check whether inputs are valid arrays, ex. GalaxiesFaintStars, dNbackground
out = super(className, self).method(args)


#########################

    for i,att in enumerate(TL.__dict__.keys()):
        print i,att

#########################

figure()
Rmask = SU.Rmask
x=((SU.Mp/const.M_earth).decompose().value);
y=((SU.Rp/const.R_earth).decompose().value);
scatter(x,y)
minorticks_on()
yscale('log');
xscale('log');
xlim(1e-1,4e5);ylim(4e-1,1e2);
yerr=[(-SU.Rperr2/const.R_earth).decompose().value,(SU.Rperr1/const.R_earth).decompose().value]
errorbar((x[Rmask]),(y[Rmask]),yerr=yerr,color='red',fmt='o')
xlabel('$M/M_\oplus$', fontsize=16)
ylabel('$R/R_\oplus$', fontsize=16)
tick_params(axis='both', which='major', labelsize=14)
savefig("scatter.png", dpi=300, transparent=True);



#########################

missionStart=60634
time = Time(float(missionStart), format='mjd', scale='tai')
self = Obs

%paste
self.r_sc
self.r_sc.to('AU')
rsc = self.r_sc
x=rsc[0]
y=rsc[1]
z=rsc[2]
c = SkyCoord(x=x,y=y,z=z, unit='AU', frame='icrs', representation='cartesian')
c
x
c.heliocentrictrueecliptic

###########################

            x = pickle.load(open(filename, 'rb'))
            if 'Name' in x:
                print x['Name']

                for attr in atts:
                    if attr in x:
                        print attr
                        
                        
                        if attr != 'radeg' or attr != 'decdeg':
                            setattr(self, attr, np.array(x[attr]))

print TL.nStars
print SU.nPlans
SS.run_sim()


DRM=SS.DRM

nobs = len(DRM)
ds = [DRM[i]['det_status'] for i in range(nobs)]
det = np.array([int(ds[i] >0) for i in range(nobs)])
ndet = np.sum(det)
cs = [DRM[i].get('char_1_time',np.nan) for i in range(nobs)]
char = np.array([int(cs[i] >0) for i in range(nobs)])
nchar = np.sum(char)

pInds = [DRM[i]['plan_inds'] for i in range(len(DRM))]
pWA = np.array([DRM[i].get('det_WA',np.nan) for i in range(len(DRM))])*u.mas


###########

figure(11)
for jitter in [0.4, 0.8, 1.6]:
    for post in [10,30]:
        pthC = '/Users/cdelacroix/Documents/AFTA-coronagraph/hlc20140623-139/hlc_20140623-139_polx_'+str(jitter)+'mas_jitter_'+str(post)+'x_contrast.fits'
        dat = fits.open(pthC)[0].data
        WA = dat[0] if dat.shape[0] == 2 else dat[:,0]
        y = dat[1] if dat.shape[0] == 2 else dat[:,1]
        plot(WA,y)
yscale('log')
figure(12)
for jitter in [0.4, 0.8, 1.6]:
    for post in [10,30]:
        pthT = '/Users/cdelacroix/Documents/AFTA-coronagraph/hlc20140623-139/hlc_20140623-139_polx_'+str(jitter)+'mas_jitter_'+str(post)+'x_throughput.fits'
        dat = fits.open(pthT)[0].data
        WA = dat[0] if dat.shape[0] == 2 else dat[:,0]
        y = dat[1] if dat.shape[0] == 2 else dat[:,1]
        plot(WA,y)

cnt=0
for jitter in [0.4, 0.8, 1.6]:
    for post in [10,30]:
        cnt+=1
        figure(cnt)
        pthP = '/Users/cdelacroix/Documents/AFTA-coronagraph/hlc20140623-139/hlc_20140623-139_polx_'+str(jitter)+'mas_jitter_'+str(post)+'x_PSF.fits'
        hdr = fits.open(pthP)[0].header
        dat = fits.open(pthP)[0].data
        imshow(dat)




