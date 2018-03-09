# -*- coding: utf-8 -*-
import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
from EXOSIMS.Prototypes.StarCatalog import StarCatalog

class KnownRVPlanetsCatalog(StarCatalog):
    """
    Star catalog based on population of known RV planets from IPAC.
    Intended for use with the KnownRVPlanets family of modules.
    
    Args: 
        \*\*specs: 
            user specified values
            
    Attributes: 
                
    
    Notes:  
    
    """

    def __init__(self, **specs):
        
                             'Spec':'st_spstr',
                             'parx':'st_plx',
                             'Umag':'st_uj',
                             'Bmag':'st_bj',
                             'Vmag':'st_vj',
                             'Rmag':'st_rc',
                             'Imag':'st_ic',
                             'Jmag':'st_j',
                             'Hmag':'st_h',
                             'Kmag':'st_k',
                             'dist':'st_dist',
                             'BV':'st_bmvj',
                             'L':'st_lum', #log(Lsun)
                             'pmra':'st_pmra', #mas/year
                             'pmdec':'st_pmdec', #mas/year
                             'rv': 'st_radv'}
        
        PPop = self.PlanetPopulation
        OS = self.OpticalSystem
        
        # IPAC data table loaded in the Planet Population module
        data = PPop.allplanetdata[:]
        
        # list of non-astropy attributes
        self.Name = 'pl_hostname',

        StarCatalog.__init__(self, ntargs=len(data), **specs)
        
        # list of astropy attributes
        self.dist = data['st_dist'].data*u.pc
        self.parx = self.dist.to('mas',equivalencies=u.parallax())
        self.coords = SkyCoord(ra=data['ra']*u.deg, dec=data['dec']*u.deg, distance=self.dist)
        self.pmra = data['st_pmra'].data*u.mas/u.yr
        self.pmdec = data['st_pmdec'].data*u.mas/u.yr


        self.rv = self.rv*u.km/u.s
        
        for att in self.atts_mapping.keys():
            ma = data[self.atts_mapping[att]]
            if type(ma.fill_value) == np.float64:
                setattr(self, att, ma.filled(np.nanmedian(ma)))
            else:
                setattr(self, att, ma.data)


    def populate_target_list(self, **specs):
        
        PPop = self.PlanetPopulation
        OS = self.OpticalSystem
        
        # filter out targets with planets outside of WA range 
        dist = data['st_dist'].filled()*u.pc
        mask = ~data['st_dist'].mask & \
                (np.arctan(PPop.sma*(1+PPop.eccen)/dist) > OS.IWA) & \
                (np.arctan(PPop.sma*(1-PPop.eccen)/dist) < OS.OWA)
        data = data[mask]
        # filter out redundant targets
        data = data[np.unique(data['pl_hostname'].data,return_index=True)[1]]
        # filter missing Vmag and BV , for OS.calc_maxintTime
        data = data[~data['st_vj'].mask]
        data = data[~data['st_bmvj'].mask]
        
        self.nStars = len(data)
        assert self.nStars, "Target list is empty: nStars = %r"%self.nStars
        

        
        self.BC =  -2.5*self.L - 26.832 - self.Vmag
        self.L = 10.**self.L
        self.MV = self.Vmag  - 5*(np.log10(self.dist.value) - 1)
        self.Binary_Cut = np.zeros(self.nStars,dtype=bool)
        
        # populate completeness values
        self.comp0 = self.Completeness.target_completeness(self)
        # populate maximum integration time
        self.maxintTime = OS.calc_maxintTime(self)
        # calculate 'true' and 'approximate' stellar masses
        self.stellar_mass()
        
        # include new attributes to the target list catalog attributes
        self.catalog_atts.append('comp0')
        self.catalog_atts.append('maxintTime')

    def filter_target_list(self, **specs):
        """ Filtering is done as part of populating the table, so this helper function
        is just a dummy."""
        
        pass
