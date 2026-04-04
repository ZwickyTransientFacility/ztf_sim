"""Code for logging observations to a sqlite database."""

import os.path
from collections import defaultdict
import uuid
import numpy as np
import pandas as pd
from sqlalchemy import create_engine, inspect, text
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import astroplan.moon
from .Fields import Fields
from .utils import *
from .constants import BASE_DIR, FILTER_ID_TO_NAME, EXPOSURE_TIME, READOUT_TIME


class ObsLogger(object):
    """Observation logger that writes pointings to a SQLite database.

    The database schema follows the 2017 LSST OpSim format with ZTF-specific
    additions. All history is also mirrored in the in-memory ``self.history``
    DataFrame for fast cadence queries.

    Attributes
    ----------
    log_name : str
        Base name of the output database file (without ``.db`` extension).
    survey_start_time : astropy.time.Time
        Reference epoch; ``expDate`` is stored as seconds elapsed since this
        time.
    history : pandas.DataFrame
        In-memory copy of the ``Summary`` table, updated after each
        `log_pointing` call.
    engine : sqlalchemy.engine.Engine
        SQLAlchemy connection to the output SQLite database.
    """

    def __init__(self, log_name, survey_start_time = Time('2018-01-01'),
            output_path = BASE_DIR+'../sims/',
            clobber = False):
        """Open (or create) the observation log database.

        Parameters
        ----------
        log_name : str
            Base name for the output SQLite file. The file is written to
            ``output_path/{log_name}.db``.
        survey_start_time : astropy.time.Time, optional
            Survey epoch used to compute ``expDate`` (seconds elapsed).
            Default is 2018-01-01.
        output_path : str, optional
            Directory in which to write the database. Default is
            ``../sims/`` relative to the package root.
        clobber : bool, optional
            If ``True``, drop and recreate the ``Field`` and ``Summary``
            tables on open. Default is ``False``.
        """
        self.log_name = log_name
        self.survey_start_time = survey_start_time
        self.prev_obs = None
        self.mjd_tonight = None
        self.moon_illumination_tonight = None
        self.engine = create_engine(
                'sqlite:///'+os.path.join(output_path,f'{self.log_name}.db'))
        self.conn = self.engine.connect()
        self.create_fields_table(clobber=clobber)
        self.create_pointing_log(clobber=clobber)

        self.history = pd.read_sql('Summary', self.engine)

    def create_fields_table(self, clobber=True):
        """Create (or recreate) the ``Field`` reference table.

        Populates the table from the ZTF field grid via `Fields`. If the
        table already exists and *clobber* is ``False``, the method does
        nothing.

        Parameters
        ----------
        clobber : bool, optional
            If ``True``, drop the existing ``Field`` table before creating a
            new one. Default is ``True``.
        """

        if clobber:
            # Drop table if it exists
            try:
                self.conn.execute("""DROP TABLE Field""")
            except:
                pass

        # If the table doesn't exist, create it
        if not inspect(self.engine).has_table('Field'): 

            self.conn.execute(text("""
            CREATE TABLE Field(
            fieldID   INTEGER PRIMARY KEY,
            fieldFov  REAL,
            fieldRA   REAL,
            fieldDec  REAL,
            fieldGL   REAL,
            fieldGB   REAL,
            fieldEL   REAL,
            fieldEB   REAL
            )"""))

            f = Fields()
            df = f.fields.reset_index()
            df.rename(columns={'field_id': 'fieldID',
                               'ra': 'fieldRA',
                               'dec': 'fieldDec',
                               'l': 'fieldGL',
                               'b': 'fieldGB',
                               'ecliptic_lon': 'fieldEL',
                               'ecliptic_lat': 'fieldEB'}, inplace=True)
            df.set_index(['fieldID'], inplace=True)
            df['fieldFov'] = 10.428

            df_min = df[['fieldFov','fieldRA', 'fieldDec', 'fieldGL', 'fieldGB',
                'fieldEL', 'fieldEB']]

            # (circumscribed) field diameter in degrees
            df_min.to_sql('Field', self.engine, if_exists='replace')

    def create_pointing_log(self, clobber=True):
        """Create (or recreate) the ``Summary`` pointing log table.

        The table schema follows the 2017 LSST OpSim ``Summary`` format with
        additional ZTF-specific columns (``totalRequestsTonight``,
        ``metricValue``, ``subprogram``).

        Parameters
        ----------
        clobber : bool, optional
            If ``True``, drop the existing ``Summary`` table before creating
            a new one. Default is ``True``.
        """

        if clobber:
            # Drop table if it exists
            try:
                self.conn.execute("""DROP TABLE Summary""")
            except:
                pass

        # If the table doesn't exist, create it
        if not inspect(self.engine).has_table('Summary'): 

            # create table
            self.conn.execute(text("""
            CREATE TABLE Summary(
            obsHistID         INTEGER PRIMARY KEY,
            requestID INTEGER,
            propID INTEGER,
            fieldID      INTEGER,
            fieldRA      REAL,
            fieldDec      REAL,
            filter             TEXT,
            expDate            INTEGER,
            expMJD             REAL,
            night              INTEGER,
            visitTime          REAL,
            visitExpTime       REAL,
            FWHMgeom           REAL,
            FWHMeff            REAL,
            airmass            REAL,
            filtSkyBright      REAL,
            lst                REAL,
            altitude           REAL,
            azimuth            REAL,
            dist2Moon          REAL,
            solarElong         REAL,
            moonRA             REAL,
            moonDec            REAL,
            moonAlt            REAL,
            moonAZ             REAL,
            moonPhase          REAL,
            sunAlt             REAL,
            sunAz              REAL,
            slewDist           REAL,
            slewTime           REAL,
            fiveSigmaDepth     REAL,
            totalRequestsTonight INTEGER,
            metricValue        REAL,
            subprogram         TEXT
            )"""))

    def log_pointing(self, state, request):
        """Record one completed observation to the database and in-memory history.

        Derives all auxiliary quantities (celestial coordinates, Moon/Sun
        positions, seeing at pointing, limiting magnitude) from *state* and
        *request*, then appends a row to ``self.history`` and writes it to
        the ``Summary`` table. Moon illumination is cached per night to avoid
        repeated expensive recomputations.

        Parameters
        ----------
        state : dict
            Telescope state dict as returned by
            ``TelescopeStateMachine.current_state_dict()`` *after* the
            exposure completes. Must include ``'current_time'`` and optionally
            ``'current_zenith_seeing'``.
        request : dict
            Observation specification returned by a queue manager. Required
            keys: ``'request_id'``, ``'target_program_id'``,
            ``'target_field_id'``, ``'target_ra'``, ``'target_dec'``,
            ``'target_filter_id'``, ``'target_exposure_time'``,
            ``'target_sky_brightness'``, ``'target_limiting_mag'``,
            ``'target_total_requests_tonight'``, ``'target_metric_value'``,
            ``'target_subprogram_name'``.
        """

        record = {}
        # don't use request_id here, but
        # let sqlite create a unique non-null key
        #record['obsHistID'] = request['request_id']
        # give request id its own column
        record['requestID'] = request['request_id']
        record['propID'] = request['target_program_id']
        record['fieldID'] = request['target_field_id']
        record['fieldRA'] = np.radians(request['target_ra'])
        record['fieldDec'] = np.radians(request['target_dec'])

        record['filter'] = FILTER_ID_TO_NAME[request['target_filter_id']]
        # times are recorded at start of exposure
        exposure_start = state['current_time'] - \
            request['target_exposure_time']
        # see note in utils.py
        exposure_start.delta_ut1_utc = 0.

        record['expDate'] = (exposure_start - self.survey_start_time).sec
        record['expMJD'] = exposure_start.mjd

        record['night'] = np.floor((exposure_start - self.survey_start_time).jd
                                   ).astype(int)
        record['visitTime'] = request[
            'target_exposure_time'].to(u.second).value
        record['visitExpTime'] = request[
            'target_exposure_time'].to(u.second).value

        # compute some values we will need
        sc = coord.SkyCoord(record['fieldRA'] * u.radian,
                            record['fieldDec'] * u.radian)
        altaz = skycoord_to_altaz(sc, exposure_start)

        if 'current_zenith_seeing' in state:
            pointing_seeing = seeing_at_pointing(state['current_zenith_seeing'].to(
            u.arcsec).value, altaz.alt.value)
            record['FWHMgeom'] = pointing_seeing
            record['FWHMeff'] = pointing_seeing

        record['airmass'] = altaz.secz.value
        record['filtSkyBright'] = request['target_sky_brightness']
        # despite the docs, it seems lst is stored as radians
        record['lst'] = np.radians(exposure_start.sidereal_time('apparent').to(
            u.hourangle).value/24.*360.)
        record['altitude'] = altaz.alt.to(u.radian).value
        record['azimuth'] = altaz.az.to(u.radian).value

        sun = coord.get_sun(exposure_start)
        sun_altaz = skycoord_to_altaz(sun, exposure_start)
        moon = coord.get_moon(exposure_start, P48_loc)
        moon_altaz = skycoord_to_altaz(moon, exposure_start)

        # WORKING AROUND a bug in sc.separation(moon)!
        #moon_sc = coord.SkyCoord(moon.ra,moon.dec)
        record['dist2Moon'] = moon.separation(sc).to(u.radian).value
        record['solarElong'] = sun.separation(sc).to(u.deg).value
        record['moonRA'] = moon.ra.to(u.radian).value
        record['moonDec'] = moon.dec.to(u.radian).value
        record['moonAlt'] = moon_altaz.alt.to(u.radian).value
        record['moonAZ'] = moon_altaz.az.to(u.radian).value

        # store tonight's mjd so that we can avoid recomputing moon
        # illumination, which profiling shows is weirdly expensive
        if np.floor(exposure_start.mjd) != self.mjd_tonight:
            self.moon_illumination_tonight = astroplan.moon.moon_illumination(
                # Don't use P48_loc to avoid astropy bug:
                # https://github.com/astropy/astroplan/pull/213
                # exposure_start, P48_loc) * 100.
                exposure_start) * 100.
            self.mjd_tonight = np.floor(exposure_start.mjd)

        record['moonPhase'] = self.moon_illumination_tonight

        record['sunAlt'] = sun_altaz.alt.to(u.radian).value
        record['sunAz'] = sun_altaz.az.to(u.radian).value
        if self.prev_obs is not None:
            sc_prev = coord.SkyCoord(self.prev_obs['fieldRA'] * u.radian,
                                     self.prev_obs['fieldDec'] * u.radian)
            record['slewDist'] = sc.separation(sc_prev).to(u.radian).value
            record['slewTime'] = (record['expDate'] -
                                  (self.prev_obs['expDate'] +
                                      self.prev_obs['visitTime']))
        record['fiveSigmaDepth'] = request['target_limiting_mag']

        # ztf_sim specific keywords!
        record['totalRequestsTonight'] = \
            request['target_total_requests_tonight']
        record['metricValue'] = request['target_metric_value']
        record['subprogram'] = request['target_subprogram_name'] 

        record_row = pd.DataFrame(record,index=[uuid.uuid1().hex])

        # append to our local history DataFrame
        # note that the index here will change when reloaded from the db
        self.history = pd.concat([self.history, record_row], axis=0, join='outer',
                                           sort=False)

        # write to the database
        record_row.to_sql('Summary', self.conn, index=False, if_exists='append')

#        # convert nan to SQL NULL. might be smarter to just replace the
#        # insertion method below with something smarter (pd.to_sql?)
#        for k,v in record.items():
#            try:
#                if np.isnan(v):
#                    record[k] = 'NULL'
#            except TypeError:
#                continue
#
#        # use placeholders to create the INSERT query
#        columns = ', '.join(list(record.keys()))
#        placeholders = '{' + '}, {'.join(list(record.keys())) + '}'
#        query = 'INSERT INTO Summary ({}) VALUES ({})'.format(
#            columns, placeholders)
#        query_filled = query.format(**record)
#        self.conn.execute(query_filled)


        # save record for next obs
        self.prev_obs = record

    def _mjd_filter_history(self, mjd_range):
        """Return a time-filtered slice of the observation history.

        Parameters
        ----------
        mjd_range : list of float or None
            ``[start_mjd, stop_mjd]`` (inclusive). If ``None``, the full
            history is returned.

        Returns
        -------
        pandas.DataFrame
            Filtered rows of ``self.history``.
        """

        if mjd_range is not None:
            assert mjd_range[0] <= mjd_range[1]
            w = ((self.history['expMJD'] >= mjd_range[0]) & 
                  (self.history['expMJD'] <= mjd_range[1])) 
            hist = self.history[w]
        else:
            hist = self.history

        return hist

    def _equivalent_obs(self, grp):
        """Convert grouped observations to equivalent standard-exposure counts.

        Accounts for both exposure time and per-observation readout overhead
        when computing the equivalent number of 30-s standard exposures.

        Parameters
        ----------
        grp : pandas.core.groupby.DataFrameGroupBy
            Grouped observations. Must contain columns ``'visitExpTime'``
            and ``'requestID'``.

        Returns
        -------
        collections.defaultdict
            Mapping from group key to equivalent exposure count (int).
        """

        total_exposure_time = grp['visitExpTime'].agg(np.sum)
        count_nobs = grp['requestID'].agg(len)

        # add readout overhead (but not slew)
        total_time = total_exposure_time + count_nobs * READOUT_TIME.to(u.second).value
        count_equivalent = np.round(total_time/(EXPOSURE_TIME + READOUT_TIME).to(u.second).value).astype(int).to_dict()

        # make this a defaultdict so we get zero values for new programs
        return defaultdict(int, count_equivalent)


    def count_equivalent_obs_by_program(self, mjd_range = None):
        """Count equivalent standard exposures grouped by program.

        Parameters
        ----------
        mjd_range : list of float or None, optional
            ``[start_mjd, stop_mjd]`` filter. Default is ``None`` (all time).

        Returns
        -------
        pandas.DataFrame
            Columns ``'program_id'`` and ``'n_obs'`` (equivalent exposure
            count).
        """


        hist = self._mjd_filter_history(mjd_range)

        grp = hist.groupby(['propID'])

        s = pd.Series(self._equivalent_obs(grp))
        s.index.name = 'program_id'
        s.name = 'n_obs'
        s = s.reset_index()
        return s

    def count_equivalent_obs_by_subprogram(self, mjd_range = None):
        """Count equivalent standard exposures grouped by program and subprogram.

        Parameters
        ----------
        mjd_range : list of float or None, optional
            ``[start_mjd, stop_mjd]`` filter. Default is ``None`` (all time).

        Returns
        -------
        pandas.DataFrame
            Columns ``'program_id'``, ``'subprogram_name'``, and ``'n_obs'``.
        """

        hist = self._mjd_filter_history(mjd_range)

        grp = hist.groupby(['propID','subprogram'])

        s = pd.Series(self._equivalent_obs(grp))
        if len(s):
            s.index.names = ['program_id','subprogram_name']
            s.name = 'n_obs'
            s = s.reset_index()
        else:
            # handle no history
            s = pd.DataFrame(columns = ['program_id','subprogram_name', 'n_obs'])
        return s

    def count_equivalent_obs_by_program_night(self, mjd_range = None):
        """Count equivalent standard exposures grouped by program and night.

        Parameters
        ----------
        mjd_range : list of float or None, optional
            ``[start_mjd, stop_mjd]`` filter. Default is ``None`` (all time).

        Returns
        -------
        pandas.DataFrame
            Columns ``'program_id'``, ``'night'``, and ``'n_obs'``.
        """

        hist = self._mjd_filter_history(mjd_range)

        grp = hist.groupby(['propID','night'])

        s = pd.Series(self._equivalent_obs(grp))
        s.index.names = ['program_id','night']
        s.name = 'n_obs'
        s = s.reset_index()
        return s

    def select_last_observed_time_by_field(self,
            field_ids = None, filter_ids = None,
            program_ids = None, subprogram_names = None,
            mjd_range = None):
        """Return the most recent observation time for each qualifying field.

        All non-``None`` filter arguments are applied with AND logic. Only
        fields that have been observed under the specified constraints are
        returned.

        Parameters
        ----------
        field_ids : set or list of int, optional
            Restrict to these field IDs.
        filter_ids : list of int, optional
            Restrict to these filter IDs (1 = g, 2 = r, 3 = i).
        program_ids : list of int, optional
            Restrict to these program IDs.
        subprogram_names : list of str, optional
            Restrict to these subprogram names.
        mjd_range : list of float or None, optional
            ``[start_mjd, stop_mjd]`` (inclusive).

        Returns
        -------
        pandas.DataFrame
            Indexed by ``fieldID`` with a single column ``'expMJD'``
            containing the maximum (most recent) MJD for each field. Fields
            with no matching observations are absent.
        """

        # start with "True" 
        w = self.history['expMJD'] > 0

        if field_ids is not None:
            w &= self.history['fieldID'].apply(lambda x: x in field_ids)

        if filter_ids is not None:
            filter_names = [FILTER_ID_TO_NAME[fi] for fi in filter_ids]
            w &= self.history['filter'].apply(lambda x: 
                    x in filter_names)

        if program_ids is not None:
            w &= self.history['propID'].apply(lambda x: 
                    x in program_ids)

        if subprogram_names is not None:
            w &= self.history['subprogram'].apply(lambda x: 
                    x in subprogram_names)

        if mjd_range is not None:
            assert mjd_range[0] <= mjd_range[1]
            w &= ((self.history['expMJD'] >= mjd_range[0]) & 
                  (self.history['expMJD'] <= mjd_range[1])) 

        # note that this only returns fields that have previously 
        # been observed under these constraints!
        return self.history.loc[
                w,['fieldID','expMJD']].groupby('fieldID').agg(np.max)

    def select_n_obs_by_field(self,
            field_ids = None, filter_ids = None,
            program_ids = None, subprogram_names = None,
            mjd_range = None):
        """Return the observation count for each qualifying field.

        All non-``None`` filter arguments are applied with AND logic. Only
        fields that have been observed at least once under the specified
        constraints are returned.

        Parameters
        ----------
        field_ids : set or list of int, optional
            Restrict to these field IDs.
        filter_ids : list of int, optional
            Restrict to these filter IDs.
        program_ids : list of int, optional
            Restrict to these program IDs.
        subprogram_names : list of str, optional
            Restrict to these subprogram names.
        mjd_range : list of float or None, optional
            ``[start_mjd, stop_mjd]`` (inclusive).

        Returns
        -------
        pandas.Series
            Indexed by ``fieldID``, named ``'n_obs'``, containing the
            observation count. Fields with no matching observations are absent.
        """

        # start with "True"
        w = self.history['expMJD'] > 0

        if field_ids is not None:
            w &= self.history['fieldID'].apply(lambda x: x in field_ids)

        if filter_ids is not None:
            filter_names = [FILTER_ID_TO_NAME[fi] for fi in filter_ids]
            w &= self.history['filter'].apply(lambda x:
                    x in filter_names)

        if program_ids is not None:
            w &= self.history['propID'].apply(lambda x:
                    x in program_ids)

        if subprogram_names is not None:
            w &= self.history['subprogram'].apply(lambda x:
                    x in subprogram_names)

        if mjd_range is not None:
            assert mjd_range[0] <= mjd_range[1]
            w &= ((self.history['expMJD'] >= mjd_range[0]) & 
                  (self.history['expMJD'] <= mjd_range[1])) 

        # note that this only returns fields that have previously 
        # been observed!   
        grp =  self.history.loc[
                w,['fieldID','expMJD']].groupby('fieldID')
        nobs = grp['expMJD'].agg(len)
        nobs.name = 'n_obs'

        return nobs

    def return_obs_history(self, time):
        """Return the observation history for the night containing *time*.

        Parameters
        ----------
        time : astropy.time.Time
            Any time within the night of interest. The night is defined as
            the 24-hour period starting at ``floor(time.mjd)``.

        Returns
        -------
        pandas.DataFrame
            Rows from ``self.history`` within the night, with columns
            ``'requestID'``, ``'propID'``, ``'fieldID'``, ``'fieldRA'``,
            ``'fieldDec'``, ``'filter'``, ``'expMJD'``, ``'visitExpTime'``,
            ``'airmass'``, ``'subprogram'``.
        """

        mjd_range = [np.floor(time.mjd), np.floor(time.mjd)+1.]
        w = ((self.history['expMJD'] >= mjd_range[0]) & 
                  (self.history['expMJD'] <= mjd_range[1])) 
        return self.history.loc[w, 
                ['requestID', 'propID', 'fieldID', 
                    'fieldRA', 'fieldDec', 'filter', 'expMJD', 'visitExpTime',
                    'airmass', 'subprogram']]

