"""Routines for working with the ZTF discrete field grid"""

import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time
from collections import defaultdict
import itertools
from .utils import *
from .SkyBrightness import SkyBrightness
from .constants import BASE_DIR, P48_loc, P48_slew_pars, PROGRAM_IDS, FILTER_IDS
from .constants import TIME_BLOCK_SIZE, MAX_AIRMASS, EXPOSURE_TIME, READOUT_TIME
from .constants import slew_time


class Fields(object):
    """Object providing convenience methods for the ZTF discrete field grid.

    Provides coordinate lookups, nightly altitude/azimuth precomputation,
    sky-limited observability estimates, and telescope-overhead time
    calculations for all fields in the grid.

    Attributes
    ----------
    fields : pandas.DataFrame
        Field properties indexed by ``field_id``. Columns: ``ra``, ``dec``,
        ``l``, ``b``, ``ecliptic_lon``, ``ecliptic_lat``, ``grid_id``.
    block_alt : pandas.DataFrame or None
        Altitude (deg) for each field (rows) at each 30-min block (columns)
        for the current night.
    block_az : pandas.DataFrame or None
        Azimuth (deg) for each field (rows) at each 30-min block (columns).
    observable_hours : pandas.Series or None
        Hours per night each field spends above the sky-limited limiting mag.
    Sky : SkyBrightness
        Sky brightness model used in observability calculations.
    """

    def __init__(self, field_filename=BASE_DIR + '../data/ZTF_Fields.txt'):
        """Load the ZTF field grid from disk.

        Parameters
        ----------
        field_filename : str, optional
            Path to the ZTF field grid text file (space-delimited, with a
            header row). Defaults to ``data/ZTF_Fields.txt`` relative to the
            package root.
        """
        self._load_fields(field_filename)
        self.loc = P48_loc
        self.current_block_night_mjd = None  # np.floor(time.mjd)
        self.current_blocks = None
        self.block_alt = None
        self.block_az = None
        self._current_observable_hours_night_mjd = None  # np.floor(time.mjd)
        self.observable_hours = None
        self.Sky = SkyBrightness()

    def _load_fields(self, field_filename):
        """Read the field grid file and populate ``self.fields``.

        Drops fields below Dec = -32° for computational speed. Assigns a
        ``grid_id`` (0–3) based on the numeric range of the field ID.

        Parameters
        ----------
        field_filename : str
            Path to the ZTF field grid text file. Expected columns (after
            skipping the header): ``field_id``, ``ra``, ``dec``, ``ebv``,
            ``l``, ``b``, ``ecliptic_lon``, ``ecliptic_lat``, ``number``.
        """

        df = pd.read_csv(field_filename,
            names=['field_id','ra','dec','ebv','l','b',
                'ecliptic_lon', 'ecliptic_lat', 'number'],
            sep='\s+',usecols=['field_id','ra','dec', 'l','b', 
                'ecliptic_lon', 'ecliptic_lat'],index_col='field_id',
            skiprows=1)


        # drop fields below dec of -32 degrees for speed
        # (grid_id = 0 has a row at -31.5)
        df = df[df['dec'] >= -32]

        # label the grid ids
        grid_id_boundaries = \
            {0: {'min':1,'max':999},
             1: {'min':1001,'max':1999},
             2: {'min':2001,'max':2999},
             3: {'min':3001,'max':3999}}

        # intialize with a bad int value
        df['grid_id'] = 99

        for grid_id, bounds in list(grid_id_boundaries.items()):
            w = (df.index >= bounds['min']) &  \
                    (df.index <= bounds['max'])
            df.loc[w,'grid_id'] = grid_id

        self.fields = df
        self.field_coords = self._field_coords()

    def _field_coords(self, cuts=None):
        """Return an astropy SkyCoord for the selected (or all) fields.

        Parameters
        ----------
        cuts : pandas.Series of bool or None, optional
            Boolean mask indexed by ``field_id``. If ``None``, all fields are
            included.

        Returns
        -------
        astropy.coordinates.SkyCoord
            ICRS coordinates for the selected fields.
        """
        if cuts is None:
            fields = self.fields
        else:
            fields = self.fields[cuts]
        return coord.SkyCoord(fields['ra'],
                              fields['dec'], frame='icrs', unit='deg')

    def compute_blocks(self, time, time_block_size=TIME_BLOCK_SIZE):
        """Pre-compute altitude and azimuth for all fields over the current night.

        Results are stored in ``self.block_alt`` and ``self.block_az``
        (DataFrames indexed by field_id, columned by block index). Also
        populates ``self.mean_observable_airmass``. Skips re-computation if
        the night has not changed since the last call.

        Parameters
        ----------
        time : astropy.time.Time
            Any time during the night of interest (used to identify the night).
        time_block_size : astropy.units.Quantity, optional
            Duration of each time block. Default is ``TIME_BLOCK_SIZE``
            (30 minutes).
        """

        # check if we've already computed for tonight:
        block_night = np.floor(time.mjd).astype(int)
        if self.current_block_night_mjd == block_night:
            return

        self.current_block_night_mjd = block_night

        blocks, times = nightly_blocks(time, time_block_size=time_block_size)
        self.current_blocks = blocks
        self.current_block_times = times

        alt_blocks = {}
        az_blocks = {}
        for bi, ti in zip(blocks, times):
            altaz = self.alt_az(ti)
            alt_blocks[bi] = altaz.alt
            az_blocks[bi] = altaz.az

        # DataFrames indexed by field_id, columns are block numbers
        self.block_alt = pd.DataFrame(alt_blocks)
        self.block_az = pd.DataFrame(az_blocks)

        block_airmass = altitude_to_airmass(self.block_alt)
        w = (block_airmass <= MAX_AIRMASS) & (block_airmass >= 1.0)
        # average airmass over the time we're above MAX_AIRMASS
        mean_observable_airmass = block_airmass[w].mean(axis=1)
        mean_observable_airmass.name = 'mean_observable_airmass'
        self.mean_observable_airmass = mean_observable_airmass

    def compute_observability(self, time, time_block_size=TIME_BLOCK_SIZE):
        """Compute the number of observable hours for each field tonight.

        Uses pre-computed block altitudes and the sky-limited limiting
        magnitude to count blocks where observations are feasible. Result
        stored in ``self.observable_hours``. Skips re-computation if the night
        has not changed since the last call.

        Parameters
        ----------
        time : astropy.time.Time
            Any time during the night of interest.
        time_block_size : astropy.units.Quantity, optional
            Duration of each time block. Default is ``TIME_BLOCK_SIZE``.
        """

        self.compute_blocks(time, time_block_size=time_block_size)

        block_night = np.floor(time.mjd).astype(int)
        if self._current_observable_hours_night_mjd == block_night:
            return

        lim_mags = {}
        # use pre-computed blocks
        for bi, ti in zip(self.current_blocks, self.current_block_times):
            df = self.fields.copy()
            df_alt = self.block_alt[bi]
            df_alt.name = 'altitude'
            df = df.join(df_alt, on='field_id')
            df_az = self.block_az[bi]
            df_az.name = 'azimuth'
            df = df.join(df_az, on='field_id')
            # for observability considerations it's sufficient to use one band
            fid=2
            df_limmag, df_sky = \
                compute_limiting_mag(df, ti, self.Sky, filter_id = fid)
            lim_mags[bi] = df_limmag

        df_lim = pd.DataFrame(lim_mags)

        observable_hours = (df_lim > 0).sum(axis=1) * \
            (TIME_BLOCK_SIZE.to(u.hour))
        observable_hours.name = 'observable_hours'
        self.observable_hours = observable_hours
        self._current_observable_hours_night_mjd = block_night

    def alt_az(self, time, cuts=None):
        """Return altitude and azimuth for all (or selected) fields at a given time.

        Parameters
        ----------
        time : astropy.time.Time
            Observation time.
        cuts : pandas.Series of bool or None, optional
            Boolean mask indexed by ``field_id``. If ``None``, all fields are
            returned. Supplying cuts is significantly slower because it
            prevents use of cached coordinates.

        Returns
        -------
        pandas.DataFrame
            DataFrame indexed by ``field_id`` with columns ``'alt'``
            (astropy Angle, degrees) and ``'az'`` (astropy Angle, degrees).
        """

        if cuts is None:
            index = self.fields.index
            fieldsAltAz = self.field_coords.transform_to(
                coord.AltAz(obstime=time, location=self.loc))
        else:
            # warning: specifying cuts makes this much slower
            index = self.fields[cuts].index
            fieldsAltAz = self._field_coords(cuts=cuts).transform_to(
                coord.AltAz(obstime=time, location=self.loc))

        return pd.DataFrame({'alt': fieldsAltAz.alt, 'az': fieldsAltAz.az},
                            index=index)

    def overhead_time(self, current_state, cuts=None):
        """Compute per-field overhead time from the current telescope position.

        Evaluates the HA, Dec, and dome axes independently and returns the
        maximum across all three plus ``READOUT_TIME``, reflecting the fact
        that readout and the slowest slew axis dominate.

        Parameters
        ----------
        current_state : dict
            Telescope state dict as returned by
            ``TelescopeStateMachine.current_state_dict()``. Required keys:
            ``'current_time'`` (astropy.time.Time), ``'current_ha'``
            (astropy Quantity), ``'current_dec'`` (astropy Quantity),
            ``'current_domeaz'`` (astropy Quantity).
        cuts : pandas.Series of bool or None, optional
            Boolean mask indexed by ``field_id`` to restrict computation to a
            subset of fields. If ``None``, all fields are evaluated.

        Returns
        -------
        overhead_df : pandas.DataFrame
            Single-column DataFrame indexed by ``field_id`` with column
            ``'overhead_time'`` (astropy.units.Quantity, seconds).
        altaz_df : pandas.DataFrame
            DataFrame indexed by ``field_id`` with columns ``'alt'`` and
            ``'az'`` at ``current_state['current_time']``.
        """

        if cuts is None:
            fields = self.fields
        else:
            fields = self.fields[cuts]

        df_altaz = self.alt_az(current_state['current_time'], cuts=cuts)
        df = fields.join(df_altaz)

        slews_by_axis = {'readout': READOUT_TIME}
        for axis in ['dome', 'dec', 'ha']:
            if axis == 'dome':
                current_coord = current_state['current_domeaz'].value
            if axis == 'ha':
                # convert to RA for ease of subtraction
                current_coord = HA_to_RA(current_state['current_ha'],
                                         current_state['current_time']).degree
            if axis == 'dec':
                current_coord = current_state['current_dec'].value
            coord = P48_slew_pars[axis]['coord']
            dangle = np.abs(df[coord] - current_coord)
            angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
            slews_by_axis[axis] = slew_time(axis, angle * u.deg)

        dfslews = pd.DataFrame(slews_by_axis, index=df.index)

        dfmax = dfslews.max(axis=1)
        dfmax = pd.DataFrame(dfmax)
        dfmax.columns = ['overhead_time']

        return dfmax, df_altaz

    def select_fields(self,
                      ra_range=None, dec_range=None,
                      l_range=None, b_range=None,
                      abs_b_range=None,
                      ecliptic_lon_range=None, ecliptic_lat_range=None,
                      grid_id=None,
                      observable_hours_range=None, 
                      field_ids = None):
        """Select a subset of fields based on sky-position criteria.

        All ``_range`` parameters accept a two-element list ``[min, max]``
        (inclusive). At most one of *b_range* and *abs_b_range* may be
        specified.

        Parameters
        ----------
        ra_range : list of float or None, optional
            RA range [min, max] in degrees.
        dec_range : list of float or None, optional
            Declination range [min, max] in degrees.
        l_range : list of float or None, optional
            Galactic longitude range [min, max] in degrees.
        b_range : list of float or None, optional
            Galactic latitude range [min, max] in degrees.
        abs_b_range : list of float or None, optional
            Absolute Galactic latitude range [min, max] in degrees.
        ecliptic_lon_range : list of float or None, optional
            Ecliptic longitude range [min, max] in degrees.
        ecliptic_lat_range : list of float or None, optional
            Ecliptic latitude range [min, max] in degrees.
        grid_id : int or None, optional
            Select only fields belonging to this grid ID (0–3).
        observable_hours_range : list of float or None, optional
            Observable hours range [min, max]. Requires
            ``compute_observability()`` to have been called first.
        field_ids : list of int or None, optional
            Restrict to this explicit list of field IDs.

        Returns
        -------
        pandas.Series of bool
            Boolean Series indexed by ``field_id``; ``True`` for fields that
            satisfy all specified criteria.

        Raises
        ------
        AssertionError
            If both *b_range* and *abs_b_range* are provided, or if
            *observable_hours_range* is specified but
            ``compute_observability()`` has not yet been called.
        """

        # start with a boolean True series:
        cuts = (self.fields['ra'] == self.fields['ra'])

        if field_ids is not None:
            in_ids = cuts.reset_index()['field_id'].apply(lambda x: x in field_ids)
            in_ids.index = cuts.index
            cuts &= in_ids

        if observable_hours_range is not None:
            # check that we've computed observable_hours
            assert(self.observable_hours is not None)
            fields = self.fields.join(self.observable_hours)
        else:
            fields = self.fields

        range_keys = ['ra', 'dec', 'l', 'b', 'ecliptic_lon', 'ecliptic_lat',
                      'observable_hours']

        assert((b_range is None) or (abs_b_range is None))

        for i, arg in enumerate([ra_range, dec_range, l_range, b_range,
                                 ecliptic_lon_range, ecliptic_lat_range,
                                 observable_hours_range]):
            if arg is not None:
                cuts = cuts & (fields[range_keys[i]] >= arg[0]) & \
                    (fields[range_keys[i]] <= arg[1])

        # easier cuts for Galactic/Extragalactic
        if abs_b_range is not None:
            cuts = cuts & (np.abs(fields['b']) >= abs_b_range[0]) & \
                (np.abs(fields['b']) <= abs_b_range[1])

        scalar_keys = ['grid_id']

        for i, arg in enumerate([grid_id]):
            if arg is not None:
                cuts = cuts & (fields[scalar_keys[i]] == arg)

        return cuts

    def select_field_ids(self, **kwargs):
        """Return the field IDs satisfying the given selection criteria.

        Thin wrapper around `select_fields` that returns indices rather than
        a boolean mask.

        Parameters
        ----------
        **kwargs
            Passed directly to `select_fields`.

        Returns
        -------
        pandas.Index
            Field IDs of the selected fields.
        """
        cuts = self.select_fields(**kwargs)
        return self.fields[cuts].index
