from __future__ import absolute_import

from builtins import object
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import astroplan.moon
from .Fields import Fields
from .utils import *
from .constants import BASE_DIR, FILTER_ID_TO_NAME


class ObsLogger(object):

    def __init__(self, log_name, survey_start_time = Time('2018-01-01'),
            clobber = False):
        self.log_name = log_name
        self.survey_start_time = survey_start_time
        self.prev_obs = None
        self.mjd_tonight = None
        self.moon_illumination_tonight = None
        self.engine = create_engine('sqlite:///{}../sims/{}.db'.format(
            BASE_DIR, self.log_name))
        self.conn = self.engine.connect()
        self.create_fields_table(clobber=clobber)
        self.create_pointing_log(clobber=clobber)

        self.history = pd.read_sql('Summary', engine)

    def create_fields_table(self, clobber=True):

        if clobber:
            # Drop table if it exists
            try:
                self.conn.execute("""DROP TABLE Field""")
            except:
                pass

        # If the table doesn't exist, create it
        if not self.engine.dialect.has_table(self.engine, 'Field'): 

            self.conn.execute("""
            CREATE TABLE Field(
            fieldID   INTEGER PRIMARY KEY,
            fieldFov  REAL,
            fieldRA   REAL,
            fieldDec  REAL,
            fieldGL   REAL,
            fieldGB   REAL,
            fieldEL   REAL,
            fieldEB   REAL
            )""")

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

        if clobber:
            # Drop table if it exists
            try:
                self.conn.execute("""DROP TABLE Summary""")
            # TODO: better error handling
            except:
                pass

        # If the table doesn't exist, create it
        if not self.engine.dialect.has_table(self.engine, 'Summary'): 

            # create table
            self.conn.execute("""
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
            )""")

    def log_pointing(self, state, request):

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

        record['filter'] = '\"' + \
            FILTER_ID_TO_NAME[request['target_filter_id']] + '\"'
        # times are recorded at start of exposure
        exposure_start = state['current_time'] - \
            request['target_exposure_time']
        # see note in utils.py
        exposure_start.delta_ut1_utc = 0.

        record['expDate'] = (exposure_start - self.survey_start_time).sec
        record['expMJD'] = exposure_start.mjd

        record['night'] = np.floor((exposure_start - self.survey_start_time).jd
                                   ).astype(np.int)
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
        record['subprogram'] = '\"' + \
            request['target_subprogram_name'] + '\"'

        record_row = pd.from_dict(record)

        # append to our local history DataFrame
        # note that the index here will change when reloaded from the db
        self.history = self.history.append(record_row)

        # write to the database
        record_row.to_sql('Summary', self.conn, index=False)

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


    def count_total_obs_by_subprogram(self):
        """Count of observations by program and subprogram.
        
        Returns a dict with keys (program_id, subprogram_name)"""

        # TODO: may need to allow for a date range to handle active months
        grp = self.history.groupby(['propID','subprogram'])

        return grp['requestID'].agg(len).to_dict()

