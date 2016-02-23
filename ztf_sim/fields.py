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
        
        dfmax = dfslews.max(axis=1)
        dfmax = pd.DataFrame(dfmax)
        dfmax.columns = ['overhead_time']
    
        return dfmax
         
            

        

def generate_test_field_grid(filename='../data/ZTF_fields.txt', 
    dbname='test_fields'):
    """Convert Eran's field grid to sqlite"""

    df = pd.read_table(filename,delimiter='\s+',skiprows=1,
            names = ['fieldid','ra','dec','extinction_b-v',
            'l','b',
            'ecliptic_lon','ecliptic_lat'],index_col=0)

    # insert label for offset grids
    grid = pd.Series(df.index >=
            1000,index=df.index,name='gridid',dtype=np.int8)

    df = df.join(grid)

    df_write_to_sqlite(df[['ra','dec','gridid']], dbname, index_label='fieldid')
