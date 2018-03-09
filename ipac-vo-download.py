#!/usr/bin/env python
#
# Usage:
#   run from shell, or
#   % python <this_file.py>
# You may need to put EXOSIMS in your $PYTHONPATH, e.g.,
#   % PYTHONPATH=/path/to/exomissionsim <this_file.py>

r"""KnownRVPlanets module unit tests

Michael Turmon, JPL, Apr. 2016
"""

import os
import json
import urllib
import urllib2
from collections import defaultdict
import numpy as np
import astropy.units as u
import astropy.constants as const


# map key names to converter-functions.  each key in the dictionary is the
# column-name within the CSV file, and each value is a function mapping the
# string (which is typically a number, but not always) to a numerical value.
# The value is translated into the right units by this function.
unit_map = dict(
    pl_hostname=str,
    pl_letter=str,
    pl_orbsmax=lambda x: float(x)*u.au,
    pl_orbsmaxerr1=lambda x: float(x)*u.au,
    pl_orbeccen=float,
    pl_orbeccenerr1=float,
    pl_bmasse=lambda x: (float(x)*const.M_earth),
    pl_bmasseerr1=lambda x: (float(x)*const.M_earth),
    pl_bmassprov=lambda x: (x == 'Msini')
    )

def load_rv_planets_json(fp):
    r"""Reads a JSON file and returns a two-level dict mapping hostnames + fields to values.

    Map to a value using a construction like, for example:
       d['Wolf 1061'][1]['pl_orbeccen']
    This is the orbital eccentricity of a Wolf 1061 planet, the second one in the catalog.
    If the value was not given in the catalog, the value in the dictionary structure
    is set up as None.
    """

    # basic data structure: a dictionary containing lists (lists of further
    # dictionaries, to be precise).
    d = defaultdict(list)
    # convert to JSON
    struct = json.load(fp)
    #print struct
    for planet in struct:
        # remap each "value" v of planet, which is a dict, through the appropriate
        # function in unit_map above -- but map empty strings to None
        row_remap = {k:(unit_map[k](v) if v else None) for (k,v) in planet.items()}
        # Append an entry to the list held within one slot of "d",
        # using the "hostname" as key.
        d[row_remap['pl_hostname']].append(row_remap)
    return d

def prettyprint(d, nhost=10):
    count = 0
    print 'Exoplanet Catalog:'
    for (key,list1) in d.items():
        if count > nhost: break
        count += 1
        print '  Host = %s' % key
        for planet in list1:
            print '    Planet = %s' % planet['pl_letter']
            for (pkey, pval) in planet.items():
                print '      %s => %s' % (pkey, str(pval))

def make_request():
    # the exoplanet API endpoint at IPAC
    # see: http://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html
    endpoint = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI"
    # parameters for the endpoint
    out_format = "json"
    # use fields='*' for everything
    # (but you need a converter for each field)
    fields = "pl_hostname,pl_letter,pl_orbsmax,pl_orbsmaxerr1,pl_orbeccen,pl_orbeccenerr1,pl_bmasse,pl_bmasseerr1,pl_bmassprov"

    params = {"table":"exoplanets",
              "format":out_format,
              "select":fields,
              "where":"pl_discmethod like 'Radial Velocity'"}

    response = urllib2.urlopen(endpoint, urllib.urlencode(params))
    return response

def main():
    response = make_request()
    d = load_rv_planets_json(response)
    prettyprint(d)

if __name__ == '__main__':
    main()
