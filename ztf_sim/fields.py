"""Routines for working with the ZTF discrete field grid"""

import numpy as np
import pandas as pd
from utils import *
import astropy.coordinates as coords
import astropy.units as u
from astropy.time import Time

class Fields:
    """Class for accessing field grid."""
    # TODO: consider using some of PTFFields.py code 
    def __init__(self,dbname='test_fields'):
        self._load_fields(dbname=dbname)
        self.loc = P48_loc
        # TODO: convert into Galactic & ecliptic coords and store
        
    def _load_fields(self, dbname='test_fields'):
        """Loads a field grid from ../data/{dbname}.db.  
        Expects fieldid, ra (deg), dec (deg) columns"""
        df = df_read_from_sqlite(dbname, index_col = 'fieldid')
        self.fields = df

    def field_coords(self):
        """Generate an astropy SkyCoord object for current fields"""
        return coords.SkyCoord(self.fields['ra'], 
            self.fields['dec'], frame='icrs', unit='deg')

    def alt_az(self, time):
        """return Altitude & Azimuth by field at a given time"""

        fieldsAltAz = self.field_coords().transform_to(
            coords.AltAz(obstime=time, location=self.loc))
        return pd.DataFrame({'alt':fieldsAltAz.alt, 'az':fieldsAltAz.az},
            index = self.fields.index) 

    def overhead_time(self, current_fieldid, time, fieldids = None):
        """Calculate overhead time in seconds from current position.
        
        fieldids is a subset of ids; if None, calculate for all ids"""
        # TODO: figure out appropriate treatment of dome at zenith

        df = self.fields.join(self.alt_az(time))

        row = df.loc[current_fieldid]

        slews_by_axis = {'readout':READOUT_TIME}
        for axis in ['HA', 'Dec', 'Dome']:
            coord = P48_slew_pars[axis]['coord']
            dangle = np.abs(df[coord] - row[coord])
            angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
            slews_by_axis[axis] = slew_time(axis, angle*u.deg)

        dfslews = pd.DataFrame(slews_by_axis,index=df.index)

        return dfslews.max(axis=1)
            
            

        


def generate_test_field_grid(dbname='test_fields'):
    """Rough code for creating a simple ZTF-scale field grid."""

    # camera view angles - TODO: update to actual numbers
    dx=np.deg2rad(7) # move along phi
    dy=np.deg2rad(7) # move along theta

    # the centers of the fields; begin with North Pole
    thetas = np.array([0.])
    phis = np.array([0.])

    for theta in np.arange(dy,np.pi,dy):
        dphi=dx/(2*np.pi*np.sin(theta+dy/2))*2*np.pi #the longest curve at theta+dy/2, should be more exact than above
        n=np.ceil(2*np.pi/dphi) # number of fields along phi at angle theta
        phi = np.arange(0,n)*(2*np.pi)/n

        thetas = np.append(thetas, theta*np.ones(np.size(phi)))
        phis = np.append(phis,phi)

    # South Pole
    thetas = np.append(thetas, np.pi)
    phis = np.append(phis, 0)

    fieldid = np.arange(len(thetas))
    ras = np.rad2deg(np.array(phis))
    #equatorial Dec=0, positive at northern part
    decs = -(np.rad2deg(np.array(thetas))-90)

    df = pd.DataFrame({'ra':ras,'dec':decs},
        index = fieldid)

    df_write_to_sqlite(df,dbname, index_label='fieldid')
