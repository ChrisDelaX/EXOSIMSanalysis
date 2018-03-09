import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

# input coordinate example
c = SkyCoord(-30, 85, 10, unit='deg,deg,pc')
ra = np.radians(c.ra)
dec = np.radians(c.dec)
print '\n input ICRS: \n ra = ', np.degrees(ra.value)*u.deg, '\n dec = ', np.degrees(dec.value)*u.deg

# Equations (3) from Leinert's paper, on page 5
eps = np.radians(23.439291) #J2000 obliquity of the ecliptic in degrees
sinbeta = np.sin(dec) * np.cos(eps) - np.cos(dec) * np.sin(eps) * np.sin(ra)
beta = np.arcsin(sinbeta)
cossollong = np.cos(ra)*np.cos(dec) / np.cos(beta)
sollong =  np.arccos(cossollong)
print '\n \n with Leinert98 equations: \n beta = ', np.degrees(beta), '\n sollong = ', np.degrees(sollong)

# python astropy conversion ra,dec -> lon,lat
c_helio = c.heliocentrictrueecliptic
print '\n \n with atropy SkyCoord: \n beta = ', c_helio.lat.value*u.deg, '\n sollong = ', c_helio.lon.value*u.deg
