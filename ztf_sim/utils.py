
import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.coordinates as coords
import astropy.units as u
from sqlalchemy import create_engine
from constants import *

def df_write_to_sqlite(df, dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df.to_sql(dbname, engine, if_exists='replace', **kwargs)

def df_read_from_sqlite(dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df = pd.read_sql(dbname, engine, **kwargs)

    return df

def _ptf_to_sqlite():
    """Convert observation history SQL dump from IPAC db to OpSim db format.

    https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions

    Fragile, as it assumes the formats below--not for general use."""

    # main dump
    #e.expid, e.prid, f.ptffield, f.objrad, f.objdecd, e.fid, e.obsmjd,
    #e.nid, e.exptime, e.airmass,
    #e.obslst, e.altitude, e.azimuth, e.moonra, e.moondec, e.moonalt, e.moonphas, \
    #e.windspeed, e.outrelhu
    df = pd.read_table('../data/opsim_dump.txt.gz',sep='|',
        names = ['obsHistID','propID','fieldID','fieldRA_deg','fieldDec_deg',
        'filter','expMJD','night','visitExpTime','airmass',
        'lst','altitude','azimuth','moonRA','moonDec','moonAlt','moonPhase',
        'wind','humidity'],index_col='obsHistID',
        skipfooter=1)

    # sky
    #e.expid, qa.fwhmsex, sdqa.metricvalue
    df_sky = pd.read_table('../data/opsim_dump_sky.txt.gz',sep='|',
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
    df_lim = pd.read_table('../data/opsim_dump_sky.txt.gz',sep='|',
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

    df_write_to_sqlite(df,'ptf')
    return df

def bin_ptf_obstimes(time_block_size = TIME_BLOCK_SIZE):
    """bin an input list of PTF exposure times (all filters,
    including H-alpha) into 
    blocks to use for weather analysis."""

    df = pd.read_table('../data/mjd.txt.gz',sep='|',
        names = ['expMJD'],
        skipfooter=1)
    t = Time(df['expMJD'],format='mjd',location=P48_loc)
    df['year'] = np.floor(t.decimalyear)
    df['block'] = block_index(t, time_block_size = TIME_BLOCK_SIZE)

    grp = df.groupby(['year','block'])
    nexps = grp.agg(len)
    nexps.rename(columns = {'expMJD':'nexps'},inplace=True)
    nexps['nexps'] = nexps['nexps'].astype(np.int8)

    df_write_to_sqlite(nexps,'weather_blocks')

def block_index(time, time_block_size = TIME_BLOCK_SIZE):
    """convert an astropy time object into a bin index for years broken up
    in time_block_size chunks."""
    from datetime import datetime


    # get the time at the start of each year
    year = np.floor(time.decimalyear)
    # this is an annoying conversion. blow up scalars:
    year = np.atleast_1d(year)
    tyear = Time([datetime(y,1,1) for y in year.astype(np.int)])

    # mjd to bin
    block_size = time_block_size.to(u.min).value
    convert = (1*u.day.to(u.min)) / block_size 

    return np.floor( (time.mjd - tyear.mjd) * convert).astype(np.int)
