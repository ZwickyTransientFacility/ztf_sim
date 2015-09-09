
import numpy as np
import pandas
from astropy.time import Time
import astropy.coordinates as coords
import astropy.units as u
from sqlalchemy import create_engine

P48_loc = coords.EarthLocation(lat=coords.Latitude('33d21m26.35s'),
    lon=coords.Longitude('-116d51m32.04s'), height=1707.)

def df_write_to_sqlite(df, dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df.to_sql(dbname, engine, if_exists='replace', **kwargs)

def df_read_from_sqlite(dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df = pandas.read_sql(dbname, engine, **kwargs)

    return df

def _ptf_to_sqlite():
    """Convert observation history SQL dump from IPAC db to OpSim db format.

    https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions

    Fragile, as it assumes the formats below--not for general use."""

    import pandas

    # main dump
    #e.expid, e.prid, f.ptffield, f.objrad, f.objdecd, e.fid, e.obsmjd,
    #e.nid, e.exptime, e.airmass,
    #e.obslst, e.altitude, e.azimuth, e.moonra, e.moondec, e.moonalt, e.moonphas, \
    #e.windspeed, e.outrelhu
    df = pandas.read_table('../data/opsim_dump.txt.gz',sep='|',
        names = ['obsHistID','propID','fieldID','fieldRA_deg','fieldDec_deg',
        'filter','expMJD','night','visitExpTime','airmass',
        'lst','altitude','azimuth','moonRA','moonDec','moonAlt','moonPhase',
        'wind','humidity'],index_col='obsHistID',
        skipfooter=1)

    # sky
    #e.expid, qa.fwhmsex, sdqa.metricvalue
    df_sky = pandas.read_table('../data/opsim_dump_sky.txt.gz',sep='|',
        names = ['obsHistID','finSeeing', 'filtSkyBrightness'],
        # have to use converters otherwise filtSkyBrightness
        # becomes object and disappears in the mean; dtypes doesn't work 
        # here becuase of skipfooter
        converters = {'filtSkyBrightness':lambda x: np.float(x)},
        skipfooter=1)

    # we have sky values and seeing on a per-CCD basis; average
    grp_sky = df_sky.groupby('obsHistID')
    seeing = grp_sky.agg(np.mean)['finSeeing']
    sky = grp_sky.agg(np.mean)['filtSkyBrightness']


    # limmag
    # e.expid, qa.fwhmsex, sdqa.metricvalue
    df_lim = pandas.read_table('../data/opsim_dump_sky.txt.gz',sep='|',
        names = ['obsHistID','finSeeing', 'fiveSigmaDepth'],
        converters = {'fiveSigmaDepth':lambda x: np.float(x)},
        skipfooter=1)

    # average by CCD
    grp_lim = df_lim.groupby('obsHistID')
    depth = grp_lim.agg(np.mean)['fiveSigmaDepth']

    df = df.join(seeing,how='outer')
    df = df.join(sky,how='outer')
    df = df.join(depth,how='outer')

    # a small number of bad rows came through
    wgood = np.isfinite(df['expMJD'])
    df = df.ix[wgood]

    # add additional columns in OpSim db
    df['sessionID'] = 0
    df['fieldRA'] = np.radians(df['fieldRA_deg'])
    df['fieldDec'] = np.radians(df['fieldDec_deg'])
    df.drop('fieldRA_deg', axis=1, inplace=True)
    df.drop('fieldDec_deg', axis=1, inplace=True)

    t = Time(df['expMJD'],format='mjd',location=P48_loc)
    df['expDate'] = t.unix - t[0].unix
    df['rotSkyPos'] = 0.
    df['ditheredRA'] = 0.
    df['ditheredDec'] = 0.

    return df

