######
import numpy as np
import pandas,os
import EXOSIMS.PlanetPhysicalModel.Forecaster
import astropy.units as u
import astropy.constants as const
import scipy.interpolate

basedir = os.path.join(os.getenv('HOME'),'INSTRUMENTS','WFIRST')
corondata = pandas.read_fwf(os.path.join(basedir,'hlc-20160125','hlc-20160125_0.4mas_jitter_1.0mas_star_results.txt'))
plandata = pandas.io.excel.read_excel(os.path.join(basedir,'WFIRST_prime_RV_targets_082016.xlsx'),header=0,skiprows=[1])

fcstr = EXOSIMS.PlanetPhysicalModel.Forecaster.Forecaster()

planradii = [fcstr.calc_radius_from_mass(np.array([v]*1000)*u.M_jupiter).mean().to('km').value for v in plandata['MSINI'].values]
planradii = np.array(planradii)*u.km

dist = plandata['DIST'].values*u.pc

#all at 65 degree phase with s = a
beta = 65*u.deg
s = plandata['A'].values*u.AU #assume s = a
r = s/np.sin(beta)
pPhib = plandata['A 65 deg'].values

fluxRatio = pPhib*((planradii/r).decompose())**2.
WA = np.arctan((s/dist).decompose()).to('arcsec')

#inside IWA - put just past IWA
iiwa = np.where(WA.value < corondata['r(arcsec)'].min())[0]
WA[iiwa] = corondata['r(arcsec)'][0:2].mean()*u.arcsec
s[iiwa] = (np.tan(WA[iiwa])*dist[iiwa]).to('AU')
r = s/np.sin(beta)
fluxRatio = pPhib*((planradii/r).decompose())**2.

#outside IWA - put just inside of OWA
oiwa = np.where(WA.value > corondata['r(arcsec)'].max())[0]
WA[oiwa] = corondata['r(arcsec)'][-2:].mean()*u.arcsec
s[oiwa] = (np.tan(WA[oiwa])*dist[oiwa]).to('AU')
r = s/np.sin(beta)
fluxRatio = pPhib*((planradii/r).decompose())**2.

#coronagraph values
speckleint = scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['I'].values)
contr = scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['contrast'].values)
corethru =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['core_thruput'].values)
psfpeak =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['PSF_peak'].values)
fwhmarea =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['area(sq_arcsec)'].values)
occtrans = psfpeak =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['occ_trans'].values)

I = speckleint(WA.to('arcsec').value)
contrast = contr(WA.to('arcsec').value)
core_thruput = corethru(WA.to('arcsec').value)
PSF_peak = psfpeak(WA.to('arcsec').value)
area = fwhmarea(WA.to('arcsec').value)
occ_trans = occtrans(WA.to('arcsec').value)


out = pandas.DataFrame({'Name':plandata['NAME'].values,
                        'Msini (M_jup)':plandata['MSINI'].values,
                        'sma (AU)':plandata['A'].values,
                        'dist (pc)':plandata['DIST'].values,
                        'p*Phi(65 deg) 565 nm':pPhib,
                        'star V mag':plandata['V'].values,
                        'prob R (R_jup)':(planradii/u.R_jup).decompose().value,
                        's (AU)':s.to('AU').value,
                        'r (AU)':r.to('AU').value,
                        'WA (as)': WA.to('arcsec').value,
                        'flux ratio':fluxRatio,
                        'I':I,
                        'contrast':contrast,
                        'core_thruput':core_thruput,
                        'PSF_peak':PSF_peak,
                        'area':area,
                        'occ_trans':occ_trans})

cols = out.columns[np.hstack((2,1,0,range(3,len(out.columns))))]
out.to_csv('RV_values_for_comparison2.txt',sep='\t',columns=cols)


######## integration time calc starts here
import pandas,os
import numpy as np
import EXOSIMS,json
import astropy.units as u
import astropy.constants as const
import scipy.interpolate


basedir = os.path.join(os.getenv('HOME'),'Documents','AFTA-coronagraph')
data = pandas.read_csv(os.path.join(basedir,'RV_values_for_comparison2.txt'),sep='\t')


#########
#corondata = pandas.read_fwf(os.path.join(basedir,'hlc-20160125','hlc-20160125_0.4mas_jitter_1.0mas_star_results.txt'))
#dist = data['dist (pc)'].values*u.pc
#planradii = data['prob R (R_jup)'].values*u.R_jupiter

#all at 65 degree phase with s = a
#beta = 65*u.deg
#s = data['sma (AU)'].values*u.AU #assume s = a
#r = s/np.sin(beta)
#pPhib = data['p*Phi(65 deg) 565 nm'].values

#fluxRatio = pPhib*((planradii/r).decompose())**2.
#WA = np.arctan((s/dist).decompose()).to('arcsec')

#outside OWA - put just inside of OWA
#oiwa = np.where(WA.value > corondata['r(arcsec)'].max())[0]
#WA[oiwa] = corondata['r(arcsec)'][-2:].mean()*u.arcsec
#s[oiwa] = (np.tan(WA[oiwa])*dist[oiwa]).to('AU')
#r = s/np.sin(beta)
#fluxRatio = pPhib*((planradii/r).decompose())**2.

#coronagraph values
#speckleint = scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['I'].values,bounds_error=False)
#contr = scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['contrast'].values,bounds_error=False)
#corethru =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['core_thruput'].values,bounds_error=False)
#psfpeak =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['PSF_peak'].values,bounds_error=False)
#fwhmarea =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['area(sq_arcsec)'].values,bounds_error=False)
#occtrans = psfpeak =  scipy.interpolate.interp1d(corondata['r(arcsec)'].values,corondata['occ_trans'].values,bounds_error=False)

#I = speckleint(WA.to('arcsec').value)
#contrast = contr(WA.to('arcsec').value)
#core_thruput = corethru(WA.to('arcsec').value)
#PSF_peak = psfpeak(WA.to('arcsec').value)
#area = fwhmarea(WA.to('arcsec').value)
#occ_trans = occtrans(WA.to('arcsec').value)


#data['s (AU)'] = s.to('AU').value
#data['r (AU)'] = r.to('AU').value
#data['WA (as)'] = WA.to('arcsec').value
#data['flux ratio'] = fluxRatio
#data['I'] = I
#data['contrast'] = contrast
#data['core_thruput'] = core_thruput
#data['PSF_peak'] = PSF_peak
#data['area'] = area
#data['occ_trans'] = occ_trans
#######


scriptfile = os.path.join(basedir,'WFIRST_KnownRV_forComparison.json')
import EXOSIMS.Prototypes.OpticalSystem
self = EXOSIMS.Prototypes.OpticalSystem.OpticalSystem(**(json.loads(open(scriptfile).read())))

inst = self.Imager
lam =  inst['lam']                             # central wavelength (given)
deltaLam = self.observingModes[0]['deltaLam']  # bandwidth (10% filter given)
QE = inst['QE'](lam)                       # quantum efficiency (ipac site)

mV = data['star V mag'].values                 # star magnitude
dMag = -2.5*np.log10(data['flux ratio'].values) # planet deltaMags

Omega = data['area'].values*u.arcsec**2     #area of FWHM in sq arcsec from PROPER
pixScale = 0.41 #given pixel scale (from Bijan's email in l/D/pix)
theta = (pixScale*lam/self.pupilDiam).decompose()*u.rad #plate scale

#!!!!!self.attenuation? 49%
#filter = 90%
#QE = 72
#!!!!!Bijan's QE doesn't match IPAC site 
#22 mag/as^2 * 1/r^2 in AU exozodi, 23 for solar sys.  one zodi each

Npix = (Omega/theta**2).decompose()   

C_F0 = self.F0(lam)*QE*self.pupilArea*deltaLam*self.attenuation #zero-mag flux (counts/s)
C_sr = (C_F0*10.**(-0.4*mV)*data['I'].values*Omega/theta**2).decompose()             # residual suppressed starlight (coro)
C_p = C_F0*10.**(-0.4*(mV + dMag))*data['core_thruput'].values         # planet signal

#C_zl = C_F0*(fZ+fEZ)*Omega                  # zodiacal light = local + exo
C_z = C_F0*10.**(-0.4*23)*Omega*data['occ_trans'].values/u.arcsec**2.
M = mV - 5*(np.log10(data['dist (pc)'].values) - 1)
C_ez = C_F0*10.**(-0.4*22)*10.**(-0.4*(M - 4.83))/(data['r (AU)'].values)**2.*Omega*data['core_thruput'].values/u.arcsec**2.
C_zl = C_z+C_ez

C_dc = Npix*inst['idark']                   # dark current
C_cc = Npix*inst['CIC']/inst['texp']        # clock-induced-charge
C_rn = Npix*(inst['sread']/inst['Gem'])**2/inst['texp'] # readout noise
C_b = inst['ENF']**2*(C_sr+C_zl+C_dc+C_cc)+C_rn     # background
ppFact = 1./30.           # post-processing contrast factor ### just parametrize by WA
C_sp = C_sr*ppFact                          # spatial structure to the speckle

SNR = 5.
intTime = (C_b+C_p)/((C_p/SNR)**2. - C_sp**2)

data = data.join(pandas.DataFrame({'intTime (hours)':intTime.to(u.hour).value}))

#data.to_csv('res3.tsv',sep='\t')





####tyler's
C_sr = 1.4e-2
C_p = 2.1e-1
C_z = 1.0e-3
C_ez =  2.3e-2
C_dc = 3.8e-3
C_cc =  7.9e-2 
C_rn = 0

C_zl = C_z+C_ez
C_b = inst['ENF']**2*(C_sr+C_zl+C_dc+C_cc)+C_rn     # background
ppFact = 1./30.           # post-processing contrast factor ### just parametrize by WA
C_sp = C_sr*ppFact                          # spatial structure to the speckle

SNR = 5.
intTime = (C_b+C_p)/((C_p/SNR)**2. - C_sp**2)

dtexp = (C_p + 2*C_b)/C_p**2.*SNR**2.
