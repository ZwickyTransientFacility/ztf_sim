
import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
import astroplan
from sqlalchemy import create_engine
from datetime import datetime
from constants import *


def df_write_to_sqlite(df, dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df.to_sql(dbname, engine, if_exists='replace', **kwargs)


def df_read_from_sqlite(dbname, **kwargs):
    engine = create_engine('sqlite:///../data/{}.db'.format(dbname))
    df = pd.read_sql(dbname, engine, **kwargs)

    return df


def HA_to_RA(ha, time):
    """convert hour angle to ra. """

    assert(time.location is not None)

    # TODO: astropy currently breaks on future dates due to IERS problems
    # issue 3275
    # hacky workaround from
    # https://groups.google.com/forum/#!msg/astropy-dev/N2Ug4RPU4DU/Gr5YNOANARUJ
    time.delta_ut1_utc = 0.
    LST = time.sidereal_time('apparent')

    ra = (LST - ha).to(u.deg)
    # wrap_angle isn't inherited
    ra.wrap_at(360. * u.deg)

    return ra


def RA_to_HA(ra, time):
    """convert ra to hour angle. """

    assert(time.location is not None)

    # TODO: astropy currently breaks on future dates due to IERS problems
    # issue 3275
    # hacky workaround from
    # https://groups.google.com/forum/#!msg/astropy-dev/N2Ug4RPU4DU/Gr5YNOANARUJ
    time.delta_ut1_utc = 0.
    LST = time.sidereal_time('apparent')

    ha = (LST - ra).to(u.deg)
    # wrap_angle isn't inherited
    ha.wrap_at(360. * u.deg)

    return ha


def previous_12deg_evening_twilight(time):
    return P48_Observer.twilight_evening_nautical(time, which='previous')


def next_12deg_evening_twilight(time):
    return P48_Observer.twilight_evening_nautical(time, which='next')


def next_12deg_morning_twilight(time):
    return P48_Observer.twilight_morning_nautical(time, which='next')


def skycoord_to_altaz(skycoord, time):
    return skycoord.transform_to(coord.AltAz(obstime=time, location=P48_loc))


def airmass_to_zenith_angle(airmass):
    return np.degrees(np.arccos(1. / airmass)) * u.deg

# cf altaz.secz.value


def airmass_to_altitude(airmass):
    return 90. * u.deg - airmass_to_zenith_angle(airmass)


def zenith_angle_to_airmass(zenith_angle):
    return 1. / np.cos(np.radians(zenith_angle))


def altitude_to_airmass(altitude):
    za = 90. - altitude  # if I make 90 a Quantity I have DataFrame troubles
    return zenith_angle_to_airmass(za)


def seeing_at_zenith(pointing_seeing, altitude):
    """Convert seeing at current pointing to zenith by multiplying by X^-3/5"""
    X = altitude_to_airmass(altitude)
    return pointing_seeing * (X**(-3. / 5.))


def seeing_at_pointing(zenith_seeing, altitude):
    """Convert zenith seeing to seeing at current altitude by multiplying by
    X^3/5"""
    X = altitude_to_airmass(altitude)
    return zenith_seeing * (X**(3. / 5.))


def approx_hours_of_darkness(time, axis=coord.Angle(23.44 * u.degree),
                             latitude=P48_loc.latitude, twilight=coord.Angle(12. * u.degree)):
    """Compute the hours of darkness (greater than t degrees twilight)

    The main approximation is a casual treatment of the time since the solstice"""
#    diff = date - pd.datetime(2000, 12, 21)
#    day = diff.total_seconds() / 24. / 3600
#    doy %= 365.25

    # TODO: actually compute the most recent solstice
    solstice = Time('2008-12-21')

    diff = (time - solstice).sec * u.second
    doy = np.floor((diff % (1 * u.year)).to(u.day).value)

    # vectorize, if needed
    # if len(np.atleast1d(doy)) == 1:
    #    doy = np.array([doy])

    m = 1. - np.tan(latitude.radian) * np.tan(axis.radian *
                                              np.cos(doy * np.pi / 182.625))
    i = np.tan(twilight.radian) / np.cos(latitude.radian)
    n = m + i
    n = np.max([0, np.min([n, 2])])
    # vectorize
    # n[n > 2] = 2
    # n[n < 0] = 0
    return 24. * u.hour * (1. - np.degrees(np.arccos(1 - n)) / 180.)


def altitude_to_fwhm(altitude, filternum):
    # values from linear fit to PTF data: in
    # notebooks/plot_sky_brightness_model.ipynb

    if filternum == 1:
        return 3.258 - 0.00925 * altitude
    elif filternum == 2:
        return 3.049 - 0.0117 * altitude
    else:
        raise NotImplementedError('FWHM not implemented for this filter')


def bin_ptf_obstimes(time_block_size=TIME_BLOCK_SIZE):
    """bin an input list of PTF exposure times (all filters,
    including H-alpha) into
    blocks to use for weather analysis."""

    df = pd.read_table('../data/mjd.txt.gz', sep='|',
                       names=['expMJD'],
                       skipfooter=1)
    t = Time(df['expMJD'], format='mjd', location=P48_loc)
    df['year'] = np.floor(t.decimalyear).astype(np.int)
    df['block'] = block_index(t, time_block_size=TIME_BLOCK_SIZE)

    grp = df.groupby(['year', 'block'])
    nexps = grp.agg(len)
    nexps.rename(columns={'expMJD': 'nexps'}, inplace=True)
    nexps['nexps'] = nexps['nexps'].astype(np.int8)

    df_write_to_sqlite(nexps, 'weather_blocks')


def block_index(time, time_block_size=TIME_BLOCK_SIZE):
    """convert an astropy time object into a bin index for years broken up
    in time_block_size chunks."""

    # get the time at the start of each year
    year = np.floor(time.decimalyear)
    # this is an annoying conversion. blow up scalars:
    year = np.atleast_1d(year)
    tyear = Time([datetime(y, 1, 1) for y in year.astype(np.int)])

    # mjd to bin
    block_size = time_block_size.to(u.min).value
    convert = (1 * u.day.to(u.min)) / block_size

    return np.floor((time.mjd - tyear.mjd) * convert).astype(np.int)


def block_index_to_time(block, time_year, where='mid',
                        time_block_size=TIME_BLOCK_SIZE):
    """Convert a block index (or array of indicies) back into astropy Time.

    block : integer or array of integers
    time_year : astropy.Time
        any time in the current year
    where : {'start', 'mid', 'end'}
        position in block to compute the time"""

    assert (where in ['start', 'mid', 'end'])

    # get the time at the start of the year
    year = np.floor(time_year.decimalyear)
    tyear = Time([datetime(np.int(year), 1, 1)])

    # this is an annoying conversion. blow up scalars:
    block = np.atleast_1d(block).astype(np.float)

    if where == 'mid':
        block += 0.5
    if where == 'end':
        block += 1
    return tyear + block * time_block_size


def nightly_blocks(time, time_block_size=TIME_BLOCK_SIZE):
    """Return block numbers and midpoint times for a given night."""

    evening_twilight = next_12deg_evening_twilight(time)
    morning_twilight = next_12deg_morning_twilight(time)

    if evening_twilight > morning_twilight:
        # the night has already started, find previous twilight
        evening_twilight = previous_12deg_evening_twilight(time)

    block_start = block_index(evening_twilight,
                              time_block_size=time_block_size)
    block_end = block_index(morning_twilight,
                            time_block_size=time_block_size)

    blocks = np.arange(block_start, block_end + 1, 1)
    times = block_index_to_time(blocks, time)

    return blocks, times


def _ptf_to_sqlite():
    """Convert observation history SQL dump from IPAC db to OpSim db format.

    https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions

    Fragile, as it assumes the formats below--not for general use."""

    # main dump
    # e.expid, e.prid, f.ptffield, f.objrad, f.objdecd, e.fid, e.obsmjd,
    # e.nid, e.exptime, e.airmass,
    # e.obslst, e.altitude, e.azimuth, e.moonra, e.moondec, e.moonalt, e.moonphas, \
    # e.windspeed, e.outrelhu
    df = pd.read_table('../data/opsim_dump.txt.gz', sep='|',
                       names=['obsHistID', 'propID', 'fieldID', 'fieldRA_deg', 'fieldDec_deg',
                              'filter', 'expMJD', 'night', 'visitExpTime', 'airmass',
                              'lst', 'altitude', 'azimuth', 'moonRA', 'moonDec', 'moonAlt', 'moonPhase',
                              'wind', 'humidity'], index_col='obsHistID',
                       skipfooter=1)

    # sky
    # e.expid, qa.fwhmsex, sdqa.metricvalue
    df_sky = pd.read_table('../data/opsim_dump_sky.txt.gz', sep='|',
                           names=['obsHistID', 'finSeeing',
                                  'filtSkyBrightness'],
                           # have to use converters otherwise filtSkyBrightness
                           # becomes object and disappears in the mean; dtypes doesn't work
                           # here becuase of skipfooter
                           converters={
                               'filtSkyBrightness': lambda x: np.float(x)},
                           skipfooter=1)

    # we have sky values and seeing on a per-CCD basis; average
    grp_sky = df_sky.groupby('obsHistID')
    seeing = grp_sky.agg(np.mean)['finSeeing']
    sky = grp_sky.agg(np.mean)['filtSkyBrightness']

    # limmag
    # e.expid, qa.fwhmsex, sdqa.metricvalue
    df_lim = pd.read_table('../data/opsim_dump_sky.txt.gz', sep='|',
                           names=['obsHistID', 'finSeeing', 'fiveSigmaDepth'],
                           converters={
                               'fiveSigmaDepth': lambda x: np.float(x)},
                           skipfooter=1)

    # average by CCD
    grp_lim = df_lim.groupby('obsHistID')
    depth = grp_lim.agg(np.mean)['fiveSigmaDepth']

    df = df.join(seeing, how='outer')
    df = df.join(sky, how='outer')
    df = df.join(depth, how='outer')

    # a small number of bad rows came through
    wgood = np.isfinite(df['expMJD'])
    df = df.ix[wgood]

    # add additional columns in OpSim db
    df['sessionID'] = 0
    df['fieldRA'] = np.radians(df['fieldRA_deg'])
    df['fieldDec'] = np.radians(df['fieldDec_deg'])
    df.drop('fieldRA_deg', axis=1, inplace=True)
    df.drop('fieldDec_deg', axis=1, inplace=True)

    t = Time(df['expMJD'], format='mjd', location=P48_loc)
    df['expDate'] = t.unix - t[0].unix
    df['rotSkyPos'] = 0.
    df['ditheredRA'] = 0.
    df['ditheredDec'] = 0.

    df_write_to_sqlite(df, 'ptf')
    return df
