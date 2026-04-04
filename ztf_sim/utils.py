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
    """Write a DataFrame to a SQLite database table.

    Parameters
    ----------
    df : pandas.DataFrame
        Data to write.
    dbname : str
        Base name of the SQLite file (without ``.db`` extension). The file
        path is constructed as ``BASE_DIR/../{directory}/{dbname}.db``.
    tablename : str or None, optional
        Table name inside the database. Defaults to *dbname* if ``None``.
    directory : str, optional
        Subdirectory relative to the package data root. Default is
        ``'data'``.
    **kwargs
        Passed to ``pandas.DataFrame.to_sql``.
    """

    if tablename is None:
        tablename = dbname
    engine = create_engine('sqlite:///{}../{}/{}.db'.format(BASE_DIR,
        directory, dbname))
    df.to_sql(tablename, engine, if_exists='replace', **kwargs)


def df_read_from_sqlite(dbname, tablename=None,
                        directory='data', **kwargs):
    """Read a table from a SQLite database into a DataFrame.

    Parameters
    ----------
    dbname : str
        Base name of the SQLite file (without ``.db`` extension).
    tablename : str or None, optional
        Table name inside the database. Defaults to *dbname* if ``None``.
    directory : str, optional
        Subdirectory relative to the package data root. Default is
        ``'data'``.
    **kwargs
        Passed to ``pandas.read_sql``.

    Returns
    -------
    pandas.DataFrame
        Contents of the requested table.
    """

    if tablename is None:
        tablename = dbname
    engine = create_engine('sqlite:///{}../{}/{}.db'.format(BASE_DIR,
        directory, dbname))
    df = pd.read_sql(tablename, engine, **kwargs)

    return df


def HA_to_RA(ha, time):
    """Convert hour angle to right ascension at the given time.

    Uses ``RA = LST - HA`` at Palomar Observatory.

    Parameters
    ----------
    ha : astropy.units.Quantity
        Hour angle in degrees (or any angular unit).
    time : astropy.time.Time
        Observation time. If ``time.location`` is ``None``, it is set to
        ``P48_loc`` in-place.

    Returns
    -------
    astropy.coordinates.Angle
        Right ascension in degrees, wrapped to [0°, 360°].
    """

    if time.location is None:
        time.location = P48_loc

    LST = time.sidereal_time('apparent')

    ra = (LST - ha).to(u.deg)
    # wrap_angle isn't inherited
    ra = ra.wrap_at(360. * u.deg)

    return ra


def RA_to_HA(ra, time):
    """Convert right ascension to hour angle at the given time.

    Uses ``HA = LST - RA`` at Palomar Observatory.

    Parameters
    ----------
    ra : astropy.units.Quantity
        Right ascension in degrees (or any angular unit).
    time : astropy.time.Time
        Observation time. If ``time.location`` is ``None``, it is set to
        ``P48_loc`` in-place.

    Returns
    -------
    astropy.coordinates.Angle
        Hour angle in degrees, wrapped to [0°, 360°].
    """

    if time.location is None:
        time.location = P48_loc

    LST = time.sidereal_time('apparent')

    ha = (LST - ra).to(u.deg)
    # wrap_angle isn't inherited
    ha = ha.wrap_at(360. * u.deg)

    return ha


def next_10deg_evening_twilight(time):
    """Return the next 10-degree evening twilight time at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −10° at sunset.
    """
    return P48_Observer.sun_set_time(time, which='next',
                                     horizon=-10*u.degree)

def next_10deg_morning_twilight(time):
    """Return the next 10-degree morning twilight time at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −10° at sunrise.
    """
    return P48_Observer.sun_rise_time(time, which='next',
                                     horizon=-10*u.degree)

def previous_12deg_evening_twilight(time):
    """Return the most recent 12-degree (nautical) evening twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the previous Sun altitude = −12° at sunset.
    """
    return P48_Observer.twilight_evening_nautical(time, which='previous')

def next_12deg_evening_twilight(time):
    """Return the next 12-degree (nautical) evening twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −12° at sunset.
    """
    return P48_Observer.twilight_evening_nautical(time, which='next')

def next_12deg_morning_twilight(time):
    """Return the next 12-degree (nautical) morning twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −12° at sunrise.
    """
    return P48_Observer.twilight_morning_nautical(time, which='next')

def next_18deg_morning_twilight(time):
    """Return the next 18-degree (astronomical) morning twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −18° at sunrise.
    """
    return P48_Observer.twilight_morning_astronomical(time, which='next')

def previous_18deg_evening_twilight(time):
    """Return the most recent 18-degree (astronomical) evening twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the previous Sun altitude = −18° at sunset.
    """
    return P48_Observer.twilight_evening_astronomical(time, which='previous')

def next_18deg_evening_twilight(time):
    """Return the next 18-degree (astronomical) evening twilight at Palomar.

    Parameters
    ----------
    time : astropy.time.Time
        Reference time.

    Returns
    -------
    astropy.time.Time
        Time of the next Sun altitude = −18° at sunset.
    """
    return P48_Observer.twilight_evening_astronomical(time, which='next')

def is_night_remaining(time):
    """Check whether dark time still remains before the next morning twilight.

    Parameters
    ----------
    time : astropy.time.Time
        Current simulation time.

    Returns
    -------
    bool
        ``True`` if *time* is before the next 12-degree morning twilight for
        the night beginning at ``floor(time.mjd)``.
    """
    Time_night_start = Time(np.floor(time.mjd), format='mjd')
    morning_twilight = next_12deg_morning_twilight(Time_night_start)

    return time <= morning_twilight


def skycoord_to_altaz(skycoord, time):
    """Transform an astropy SkyCoord to altitude–azimuth at Palomar.

    Parameters
    ----------
    skycoord : astropy.coordinates.SkyCoord
        ICRS sky coordinate.
    time : astropy.time.Time
        Observation time.

    Returns
    -------
    astropy.coordinates.SkyCoord
        Coordinate in the AltAz frame at ``P48_loc``.
    """
    return skycoord.transform_to(coord.AltAz(obstime=time, location=P48_loc))


def airmass_to_zenith_angle(airmass):
    """Convert airmass to zenith angle.

    Uses the plane-parallel approximation ``X = 1 / cos(z)``.

    Parameters
    ----------
    airmass : float or array-like
        Airmass value(s).

    Returns
    -------
    astropy.units.Quantity
        Zenith angle in degrees.
    """
    return np.degrees(np.arccos(1. / airmass)) * u.deg

# cf altaz.secz.value


def airmass_to_altitude(airmass):
    """Convert airmass to altitude above the horizon.

    Parameters
    ----------
    airmass : float or array-like
        Airmass value(s).

    Returns
    -------
    astropy.units.Quantity
        Altitude in degrees (= 90° − zenith angle).
    """
    return 90. * u.deg - airmass_to_zenith_angle(airmass)


def zenith_angle_to_airmass(zenith_angle):
    """Convert zenith angle to airmass using the plane-parallel approximation.

    Parameters
    ----------
    zenith_angle : float or array-like
        Zenith angle in degrees.

    Returns
    -------
    float or numpy.ndarray
        Airmass ``X = 1 / cos(zenith_angle)``.
    """
    return 1. / np.cos(np.radians(zenith_angle))


def altitude_to_airmass(altitude):
    """Convert altitude to airmass using the plane-parallel approximation.

    Parameters
    ----------
    altitude : float or array-like
        Pointing altitude in degrees.

    Returns
    -------
    float or numpy.ndarray
        Airmass ``X = 1 / cos(90° − altitude)``.
    """
    za = 90. - altitude  # if I make 90 a Quantity I have DataFrame troubles
    return zenith_angle_to_airmass(za)

def maximum_altitude(dec, lat=P48_loc.lat.degree):
    """Compute the transit altitude for a source at the given declination.

    Parameters
    ----------
    dec : float or pandas.Series
        Source declination in degrees.
    lat : float, optional
        Observer latitude in degrees. Default is the P48 latitude.

    Returns
    -------
    float or pandas.Series
        Transit altitude in degrees. Same type and shape as *dec*.

    Notes
    -----
    Let ``px = 90 - dec`` and ``pz = 90 - lat``.

    * If ``px >= pz`` (source transits south of zenith):
      ``alt = 90 - lat + dec``
    * If ``px < pz`` (source transits north of zenith):
      ``alt = 90 + lat - dec``
    """

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
    """Convert seeing at the pointing altitude to zenith seeing.

    Parameters
    ----------
    pointing_seeing : float or array-like
        Seeing FWHM at the pointing altitude (arcsec).
    altitude : float or array-like
        Pointing altitude in degrees.

    Returns
    -------
    float or numpy.ndarray
        Zenith seeing FWHM in the same units as *pointing_seeing*.
    """
    X = altitude_to_airmass(altitude)
    return pointing_seeing * (X**(-3. / 5.))


def seeing_at_pointing(zenith_seeing, altitude):
    """Convert zenith seeing to seeing at the pointing altitude.

    Parameters
    ----------
    zenith_seeing : float or array-like
        Zenith seeing FWHM (arcsec, or any consistent unit).
    altitude : float or array-like
        Pointing altitude in degrees.

    Returns
    -------
    float or numpy.ndarray
        Seeing FWHM at the pointing altitude, in the same units as
        *zenith_seeing*.
    """
    X = altitude_to_airmass(altitude)
    return zenith_seeing * (X**(3. / 5.))


def approx_hours_of_darkness(time, axis=coord.Angle(23.44 * u.degree),
                             latitude=P48_loc.lat, twilight=coord.Angle(12. * u.degree)):
    """Estimate the hours of darkness (Sun below −twilight) for a given night.

    Parameters
    ----------
    time : astropy.time.Time
        Any time within the night of interest.
    axis : astropy.coordinates.Angle, optional
        Obliquity of the ecliptic (Earth's axial tilt). Default is 23.44°.
    latitude : astropy.coordinates.Angle, optional
        Observer latitude. Default is the P48 latitude.
    twilight : astropy.coordinates.Angle, optional
        Twilight depression angle. Default is 12° (nautical twilight).

    Returns
    -------
    astropy.units.Quantity
        Approximate hours of darkness as an array (same shape as *time*
        after ``np.atleast_1d``).

    Notes
    -----
    Uses a day-of-year approximation relative to the 2016 winter solstice.
    Accurate to roughly ±15 minutes for typical use; prefer
    ``astroplan`` twilight functions for precision work.
    """

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
    """Estimate seeing FWHM from pointing altitude using an empirical fit to PTF data.

    Parameters
    ----------
    altitude : float
        Pointing altitude in degrees.
    filternum : int
        Filter identifier (1 = g, 2 = r, 3 = i). g and i use the same fit
        as they have sparse PTF data.

    Returns
    -------
    float
        Predicted FWHM in arcseconds.

    Raises
    ------
    NotImplementedError
        If *filternum* is not 1, 2, or 3.

    Notes
    -----
    Coefficients from a linear fit to PTF DIQ data (see
    ``notebooks/plot_sky_brightness_model.ipynb``):

    * g / i : ``FWHM = 3.258 - 0.00925 * altitude``
    * r     : ``FWHM = 3.049 - 0.0117 * altitude``
    """
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
    """Bin PTF exposure times into blocks and write to the weather SQLite database.

    Reads PTF MJD timestamps from ``data/mjd.txt.gz``, assigns each exposure
    to a ``(year, block_index)`` bin, counts the number of exposures per bin,
    and writes the result to ``data/weather_blocks.db`` via
    `df_write_to_sqlite`. This is a one-time utility used to build the
    historical weather database consumed by `PTFObservabilityDB`.

    Parameters
    ----------
    time_block_size : astropy.units.Quantity, optional
        Duration of each time block. Default is ``TIME_BLOCK_SIZE`` (30 min).
    """

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
    """Convert an astropy Time to an integer block index within the current year.

    The block index counts 30-minute (or *time_block_size*) bins since the
    start of the calendar year.

    Parameters
    ----------
    time : astropy.time.Time
        Input time(s).
    time_block_size : astropy.units.Quantity, optional
        Block duration. Default is ``TIME_BLOCK_SIZE`` (30 min).

    Returns
    -------
    numpy.ndarray of int
        Block index for each element of *time*.

    Notes
    -----
    Formula::

        block = floor((time.mjd - year_start_mjd) * (1440 / block_size_min))
    """

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
    """Convert a block index (or array of indices) back to astropy Time.

    Parameters
    ----------
    block : int or array-like of int
        Block index within the year.
    time_year : astropy.time.Time
        Any time in the target year (only the year component is used).
    where : {'start', 'mid', 'end'}, optional
        Which part of the block to return. Default is ``'mid'``.
    time_block_size : astropy.units.Quantity, optional
        Block duration. Default is ``TIME_BLOCK_SIZE``.

    Returns
    -------
    astropy.time.Time
        Time at the requested position within the block(s).

    Raises
    ------
    AssertionError
        If *where* is not one of ``'start'``, ``'mid'``, ``'end'``.
    """

    assert (where in ['start', 'mid', 'end'])

    # get the time at the start of the year
    year = np.floor(time_year.decimalyear)
    tyear = Time([datetime(int(year), 1, 1)])

    # this is an annoying conversion. blow up scalars:
    block = np.atleast_1d(block).astype(float)

    if where == 'mid':
        block += 0.5
    if where == 'end':
        block += 1
    return tyear + block * time_block_size


def nightly_blocks(time, time_block_size=TIME_BLOCK_SIZE):
    """Return block indices and midpoint times for the night containing *time*.

    The night is defined from 12-degree evening twilight to 12-degree morning
    twilight at Palomar. If the night has already started (evening twilight is
    in the past), the previous evening twilight is used.

    Parameters
    ----------
    time : astropy.time.Time
        Any time associated with the night of interest.
    time_block_size : astropy.units.Quantity, optional
        Block duration. Default is ``TIME_BLOCK_SIZE``.

    Returns
    -------
    blocks : numpy.ndarray of int
        Block indices from evening to morning twilight (inclusive).
    times : astropy.time.Time
        Midpoint times of each block.
    """

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
    """Compute the fraction of a block covered by an observation window.

    Handles four cases: the window completely covers the block, the window
    is entirely within the block, the window starts inside and ends later,
    and the window starts earlier and ends inside.

    Parameters
    ----------
    block_index : int
        Index of the block to evaluate.
    obs_start_time : astropy.time.Time
        Start of the observation window.
    obs_end_time : astropy.time.Time
        End of the observation window.

    Returns
    -------
    float
        Fraction of the block (0–1) occupied by the observation window.

    Raises
    ------
    AssertionError
        If no case matches (should never occur for valid inputs).

    Notes
    -----
    Only scalar block indices are supported.
    """

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
    """Return ``len(np.atleast_1d(x))`` to handle scalar or array inputs.

    Parameters
    ----------
    x : scalar or array-like
        Input value.

    Returns
    -------
    int
        Length of *x* after promoting it to a 1-D array.
    """
    return len(np.atleast_1d(x))


def compute_limiting_mag(df, time, sky, filter_id=None):
    """Compute the 5σ limiting magnitude for all fields in a DataFrame.

    Predicts sky brightness with *sky*, computes the seeing at each pointing
    altitude, and calls `limiting_mag`. Applies a per-filter renormalisation
    so that g- and i-band magnitudes span the same dynamic range as r-band
    for use in the Gurobi metric. Also locks out fields that are too low
    (below MAX_AIRMASS altitude), within 20° of the Moon, or violate the
    P48 Reed pointing limits.

    Parameters
    ----------
    df : pandas.DataFrame
        Fields to evaluate. Required columns: ``ra``, ``dec``, ``altitude``,
        ``azimuth``, ``filter_id`` (unless overridden by *filter_id*).
    time : astropy.time.Time
        Block midpoint time (used for Sun/Moon position).
    sky : SkyBrightness or FakeSkyBrightness
        Sky brightness model.
    filter_id : int or None, optional
        If provided, override the ``'filter_id'`` column in *df*. Default is
        ``None`` (use the column as-is).

    Returns
    -------
    limiting_mag : pandas.Series of float
        5σ limiting AB magnitude per field. Fields violating any constraint
        receive a value of −99.
    sky_brightness : pandas.Series of float
        Predicted sky surface brightness in AB mag arcsec⁻².

    Raises
    ------
    ValueError
        If any pointing has a Sun altitude above −6° (inside 6-degree
        twilight).
    """

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
    """Convert PTF observation history from an IPAC SQL dump to OpSim format.

    Reads the main PTF dump (``data/opsim_dump.txt.gz``), sky-brightness
    dump (``data/opsim_dump_sky.txt.gz``), and limiting-magnitude dump, then
    merges them, adds OpSim-compatible derived columns, and writes the result
    to ``data/ptf.db`` via `df_write_to_sqlite`.

    This is a one-time utility function and is not part of the normal
    simulation workflow. It is fragile and assumes specific column formats.

    Returns
    -------
    pandas.DataFrame
        Merged PTF observation table.

    Notes
    -----
    Schema reference:
    https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions
    """

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
                               'filtSkyBrightness': lambda x: float(x)},
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
                               'fiveSigmaDepth': lambda x: float(x)},
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
    """Export pointing data to a space-delimited text file for the image simulator.

    Reads the ``Summary`` table from ``../sims/{dbname}.db`` and writes
    selected columns (``ra``, ``dec``, ``fieldID``, ``filter``,
    ``imagetype``, ``expMJD``) to ``../sims/{dbname}.txt``.

    Parameters
    ----------
    dbname : str
        Base name of the simulation database (without ``.db`` extension).
    **kwargs
        Passed to ``pandas.read_sql`` (e.g., ``chunksize``).
    """

    engine = create_engine('sqlite:///../sims/{}.db'.format(dbname))
    df = pd.read_sql('Summary', engine, **kwargs)

    df['ra'] = np.degrees(df['fieldRA'])
    df['dec'] = np.degrees(df['fieldDec'])
    df['imagetype'] = 0

    df[['ra', 'dec', 'fieldID',
        'filter', 'imagetype', 'expMJD']].to_csv('../sims/{}.txt'.format(
            dbname), sep=' ', header=False, index=False)
