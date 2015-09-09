"""Routines for working with the ZTF discrete field grid"""

import numpy as np
import pandas
from utils import df_write_to_sqlite, df_read_from_sqlite

def load_fields(dbname='test_fields'):
    df = df_read_from_sqlite(dbname, index_col = 'fieldid')
    return df


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

    df = pandas.DataFrame({'ra':ras,'dec':decs},
        index = fieldid)

    df_write_to_sqlite(df,dbname, index_label='fieldid')

