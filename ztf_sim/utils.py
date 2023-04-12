"""Utility routines."""

import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
import astroplan
from sqlalchemy import create_engine
from datetime import datetime
from .constants import BASE_DIR, P48_loc, P48_Observer, TIME_BLOCK_SIZE
from .constants import EXPOSURE_TIME, MAX_AIRMASS 
from .magnitudes import limiting_mag





def df_write_to_sqlite(df, dbname, tablename=None,
                       directory='data', **kwargs):

    if tablename is None:
        tablename = dbname
    engine = create_engine('sqlite:///{}../{}/{}.db'.format(BASE_DIR,
        directory, dbname))
    df.to_sql(tablename, engine, if_exists='replace', **kwargs)


def df_read_from_sqlite(dbname, tablename=None,
                        directory='data', **kwargs):

    if tablename is None:
        tablename = dbname
    engine = create_engine('sqlite:///{}../{}/{}.db'.format(BASE_DIR,
        directory, dbname))
    df = pd.read_sql(tablename, engine, **kwargs)

    return df


def HA_to_RA(ha, time):
    """convert hour angle to ra. """

    if time.location is None:
        time.location = P48_loc

    LST = time.sidereal_time('apparent')

    ra = (LST - ha).to(u.deg)
    # wrap_angle isn't inherited
    ra = ra.wrap_at(360. * u.deg)

    return ra


def RA_to_HA(ra, time):
    """convert ra to hour angle. """

    if time.location is None:
        time.location = P48_loc

    LST = time.sidereal_time('apparent')

    ha = (LST - ra).to(u.deg)
    # wrap_angle isn't inherited
    ha = ha.wrap_at(360. * u.deg)

    return ha


def previous_12deg_evening_twilight(time):
    return P48_Observer.twilight_evening_nautical(time, which='previous')


def next_12deg_evening_twilight(time):
    return P48_Observer.twilight_evening_nautical(time, which='next')

def next_12deg_morning_twilight(time):
    return P48_Observer.twilight_morning_nautical(time, which='next')

def next_18deg_morning_twilight(time):
    return P48_Observer.twilight_morning_astronomical(time, which='next')

def previous_18deg_evening_twilight(time):
    return P48_Observer.twilight_evening_astronomical(time, which='previous')

def next_18deg_evening_twilight(time):
    return P48_Observer.twilight_evening_astronomical(time, which='next')


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

def maximum_altitude(dec, lat=P48_loc.lat.degree):
    """Compute the altitude of a source with declination dec as it transits the
    meridian.
    
    dec: Pandas DataFrame, which may be Multi-Indexed"""

    px = 90 - dec
    pz = 90 - lat

    results = dec*0.

    w = (px >= pz)

    if np.sum(w.values.flatten()):
        results[w] = 90 - lat + dec[w]
    if np.sum(~w.values.flatten()):
        results[~w] = 90 + lat - dec[~w]

    return results



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
                             latitude=P48_loc.lat, twilight=coord.Angle(12. * u.degree)):
    """Compute the hours of darkness (greater than t degrees twilight)

    The main approximation is a casual treatment of the time since the solstice"""

    # would be better to actually compute a recent solstice
    solstice = Time('2016-12-21')

    diff = (time - solstice).sec * u.second
    doy = np.floor((diff % (1 * u.year)).to(u.day).value)

    # vectorize, if needed
    if len(np.atleast_1d(doy)) == 1:
        doy = np.array([doy])

    m = 1. - np.tan(latitude.radian) * np.tan(axis.radian *
                                              np.cos(doy * np.pi / 182.625))
    i = np.tan(twilight.radian) / np.cos(latitude.radian)
    n = m + i
    #n = np.max([0, np.min([n, 2])])
    # vectorize
    n[n > 2] = 2
    n[n < 0] = 0
    return 24. * u.hour * (1. - np.degrees(np.arccos(1 - n)) / 180.)


def altitude_to_fwhm(altitude, filternum):
    # values from linear fit to PTF data: in
    # notebooks/plot_sky_brightness_model.ipynb

    # don't have a lot of PTF i-band data, so let's make it the same as
    # r-band (atmosphere should contribute less)
    if (filternum == 1) or (filternum == 3):
        return 3.258 - 0.00925 * altitude
    elif filternum == 2:
        return 3.049 - 0.0117 * altitude
    else:
        raise NotImplementedError('FWHM not implemented for this filter')


def bin_ptf_obstimes(time_block_size=TIME_BLOCK_SIZE):
    """bin an input list of PTF exposure times (all filters,
    including H-alpha) into
    blocks to use for weather analysis."""

    df = pd.read_table(BASE_DIR + '../data/mjd.txt.gz', sep='|',
                       names=['expMJD'],
                       skipfooter=1)
    t = Time(df['expMJD'], format='mjd', location=P48_loc)
    df['year'] = np.floor(t.decimalyear).astype(int)
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
    tyear = Time([datetime(y, 1, 1) for y in year.astype(int)])

    # mjd to bin
    block_size = time_block_size.to(u.min).value
    convert = (1 * u.day.to(u.min)) / block_size

    return np.floor((time.mjd - tyear.mjd) * convert).astype(int)


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
    tyear = Time([datetime(int(year), 1, 1)])

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

    if ((evening_twilight > morning_twilight) 
        or (evening_twilight.value < 0)
        # some versions of astroplan+astropy seem to return masked arrays
        # instead
        or (type(evening_twilight.value) == np.ma.core.MaskedArray)):
        # the night has already started, find previous twilight
        evening_twilight = previous_12deg_evening_twilight(time)

    block_start = block_index(evening_twilight,
                              time_block_size=time_block_size)
    block_end = block_index(morning_twilight,
                            time_block_size=time_block_size)

    blocks = np.arange(block_start, block_end + 1, 1)
    times = block_index_to_time(blocks, time, where='mid')

    return blocks, times

def block_use_fraction(block_index, obs_start_time, obs_end_time):
    """Given a block index and Times specifying the start and end of an observation window, return the fraction of the block covered by the window.
    
    Scalars only for now"""

    # obs_start_time is just providing the year here
    block_tstart = block_index_to_time(block_index, obs_start_time,
            where='start')[0]
    block_tend = block_index_to_time(block_index, obs_end_time,
            where='end')[0]

    # block completely filled
    if (obs_start_time <= block_tstart) and (obs_end_time >= block_tend):
        return 1.0

    # window completely within the block
    if (obs_start_time >= block_tstart) and (obs_end_time <= block_tend):
        return ((obs_end_time - obs_start_time) /
                TIME_BLOCK_SIZE).to(u.dimensionless_unscaled).value

    # window starts within the block and finishes in a later block
    if (obs_start_time > block_tstart) and (obs_end_time > block_tend):
        return ((block_tend - obs_start_time) /
                TIME_BLOCK_SIZE).to(u.dimensionless_unscaled).value

    # window starts in an earlier block and finishes in this block
    if (obs_start_time < block_tstart) and (obs_end_time < block_tend):
        return ((obs_end_time - block_tstart) /
                TIME_BLOCK_SIZE).to(u.dimensionless_unscaled).value

    # this should never be reached
    raise AssertionError('Block use calculation is inconsistent')





def scalar_len(x):
    """Convenience function to sanitize potential scalars or arrays
    so they can return a len"""
    return len(np.atleast_1d(x))


def compute_limiting_mag(df, time, sky, filter_id=None):
    """compute limiting magnitude based on sky brightness and seeing
    
    df is a DataFrame of fields"""

    # copy df so we can edit the filter id if desired
    if filter_id is not None:
        df = df.copy()
        df['filter_id'] = filter_id

    # compute inputs for sky brightness
    sc = coord.SkyCoord(df['ra'], df['dec'], frame='icrs', unit='deg')
    sun = coord.get_sun(time)
    sun_altaz = skycoord_to_altaz(sun, time)
    moon = coord.get_moon(time, location=P48_loc)
    moon_altaz = skycoord_to_altaz(moon, time)
    df.loc[:, 'moonillf'] = astroplan.moon.moon_illumination(time)
    
    # WORKING AROUND BUG in moon distance!!!!  171110
    df.loc[:, 'moon_dist'] = moon.separation(sc).to(u.deg).value
    df.loc[:, 'moonalt'] = moon_altaz.alt.to(u.deg).value
    df.loc[:, 'sunalt'] = sun_altaz.alt.to(u.deg).value

    # check if the sun is up anywhere and break things if it isn't
    if np.sum(df['sunalt'] > -6) != 0:
        raise ValueError('Some pointings outside six-degree twilight!')

    # compute sky brightness
    # only have values for reasonable altitudes (set by R20_absorbed...)
    wup = df['altitude'] >= airmass_to_altitude(MAX_AIRMASS) 
    df.loc[wup, 'sky_brightness'] = sky.predict(df[wup])

    # compute seeing at each pointing
    df.loc[wup, 'seeing'] = seeing_at_pointing(2.0*u.arcsec, 
        df.loc[wup,'altitude'])

    df.loc[wup, 'limiting_mag'] = limiting_mag(EXPOSURE_TIME, 
        df.loc[wup, 'seeing'],
        df.loc[wup, 'sky_brightness'],
        filter_id = df.loc[wup,'filter_id'],
        altitude = df.loc[wup,'altitude'], SNR=5.)

    # renormalize limiting mags to the R-band range so we maintain 
    # the causal structure with airmass, etc. but can get i-band scheduled
    
    # bright time limiting mags (from PTF-trained model--see 170930 notes
    # and plot_sky_brightness_model.ipynb)
    mlim_bright_g = 19.9
    mlim_bright_r = 20.1
    mlim_bright_i = 19.5
    dm_g = (21.9-19.9)
    dm_r = (21.5-20.1)
    dm_i = (20.9-19.5)

    wg = df['filter_id'] == 1
    if np.sum(wg):
        df.loc[wg,'limiting_mag'] = \
            (df.loc[wg,'limiting_mag'] - mlim_bright_g) * dm_r/dm_g \
            + mlim_bright_r

    wi = df['filter_id'] == 3
    if np.sum(wi):
        df.loc[wi,'limiting_mag'] = \
            (df.loc[wi,'limiting_mag'] - mlim_bright_i) * dm_r/dm_i \
            + mlim_bright_r

    # assign a very bright limiting mag to the fields that are down 
    # so the metric goes to zero
    df.loc[~wup, 'limiting_mag'] = -99

    # assign a very bright limiting mag to the fields within 20 degrees of
    # the moon 
    wmoon = df['moon_dist'] < 20
    df.loc[wmoon, 'limiting_mag'] = -99

    # need to check the Hour Angle at both the start and the end of the
    # block, since we don't know the exact time it will be observed

    # time is provided at the block midpoint

    ha_vals = RA_to_HA(df['ra'].values*u.degree, 
            time - TIME_BLOCK_SIZE/2.)
    # for limits below, need ha-180-180
    ha_vals = ha_vals.wrap_at(180.*u.degree)
    ha = pd.Series(ha_vals.to(u.degree).value, index=df.index, name='ha')

    ha_vals_end = RA_to_HA(df['ra'].values*u.degree, 
            time + TIME_BLOCK_SIZE/2.)
    # for limits below, need ha-180-180
    ha_vals_end = ha_vals_end.wrap_at(180.*u.degree)
    ha_end = pd.Series(ha_vals_end.to(u.degree).value, index=df.index, name='ha')

    # lock out TCS limits
    
    # Reed limits |HA| to < 5.95 hours (most relevant for circumpolar
    # fields not hit by the airmass cut)
    whalimit = np.abs(ha) >= (5.95 * u.hourangle).to(u.degree).value
    whalimit_end = np.abs(ha_end) >= (5.95 * u.hourangle).to(u.degree).value
    df.loc[whalimit | whalimit_end, 'limiting_mag'] = -99

    # 1) HA < -17.6 deg && Dec < -22 deg is rejected for both track & stow because of interference with FFI.
    w1 = (ha <= -17.6) & (df['dec'] <= -22)
    w1_end = (ha_end <= -17.6) & (df['dec'] <= -22)
    df.loc[w1 | w1_end, 'limiting_mag'] = -99

    # West of HA -17.6 deg, Dec < -45 deg is rejected for tracking because of the service platform in the south.  
    w2 = (ha >= -17.6) & (df['dec'] <= -45)
    w2_end = (ha_end >= -17.6) & (df['dec'] <= -45)
    df.loc[w2 | w2_end, 'limiting_mag'] = -99

    # fabs(HA) > 3 deg is rejected for Dec < -46 to protect the shutter "ears".  
    w3 = (np.abs(ha) >= 3.) & (df['dec'] <= -46)
    w3_end = (np.abs(ha_end) >= 3.) & (df['dec'] <= -46)
    df.loc[w3 | w3_end, 'limiting_mag'] = -99

    # dec > 87.5 is rejected
    w4 = (df['dec'] > 87.5)
    df.loc[w4, 'limiting_mag'] = -99

    return df['limiting_mag'], df['sky_brightness']

def _ptf_to_sqlite():
    """Convert observation history SQL dump from IPAC db to OpSim db format.

    https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions

    Fragile, as it assumes the formats below--not for general use."""

    # main dump
    # e.expid, e.prid, f.ptffield, f.objrad, f.objdecd, e.fid, e.obsmjd,
    # e.nid, e.exptime, e.airmass,
    # e.obslst, e.altitude, e.azimuth, e.moonra, e.moondec, e.moonalt, e.moonphas, \
    # e.windspeed, e.outrelhu
    df = pd.read_table(BASE_DIR + '../data/opsim_dump.txt.gz', sep='|',
                       names=['obsHistID', 'propID', 'fieldID', 'fieldRA_deg', 'fieldDec_deg',
                              'filter', 'expMJD', 'night', 'visitExpTime', 'airmass',
                              'lst', 'altitude', 'azimuth', 'moonRA', 'moonDec', 'moonAlt', 'moonPhase',
                              'wind', 'humidity'], index_col='obsHistID',
                       skipfooter=1)

    # sky
    # e.expid, qa.fwhmsex, sdqa.metricvalue
    df_sky = pd.read_table(BASE_DIR + '../data/opsim_dump_sky.txt.gz', sep='|',
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
    df_lim = pd.read_table(BASE_DIR + '../data/opsim_dump_sky.txt.gz', sep='|',
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
    df = df.loc[wgood]

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

    df['filter'] = df['filter'].map({1:'g',2:'r',4:'i'})
    df.sort_values('expMJD',inplace=True)

    # for some reason the night values from the db are not monotonic in MJD
    # make my own versions
    df['night'] = np.floor(df['expMJD'] - 54847).astype(int)

    df_write_to_sqlite(df, 'ptf', tablename='Summary')
    return df


def export_pointings_to_surace(dbname, **kwargs):
    """put pointing data in format Jason Surace wants for his image simulator"""

    engine = create_engine('sqlite:///../sims/{}.db'.format(dbname))
    df = pd.read_sql('Summary', engine, **kwargs)

    df['ra'] = np.degrees(df['fieldRA'])
    df['dec'] = np.degrees(df['fieldDec'])
    df['imagetype'] = 0

    df[['ra', 'dec', 'fieldID',
        'filter', 'imagetype', 'expMJD']].to_csv('../sims/{}.txt'.format(
            dbname), sep=' ', header=False, index=False)
