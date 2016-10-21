
import numpy as np
from sqlalchemy import create_engine
import sqlite3
import astropy.coordinates as coord
import astropy.units as u
import astroplan.moon
from utils import *
from constants import *


class ObsLogger:

    def __init__(self, run_name, survey_start_time):
        self.run_name = run_name
        self.survey_start_time = survey_start_time
        self.prev_obs = None
        self.mjd_tonight = None
        self.moon_illumination_tonight = None
        self.engine = create_engine('sqlite:///../sims/{}.db'.format(
            self.run_name))
        self.conn = self.engine.connect()

    def create_pointing_log(self, clobber=True):

        if clobber:
            try:
                self.conn.execute("""DROP TABLE Summary""")
            # TODO: better error handling
            except:
                pass

        self.conn.execute("""
        CREATE TABLE Summary(
        obsHistID         INTEGER PRIMARY KEY,
        sessionID INTEGER,
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
        finRank            REAL,
        FWHMgeom           REAL,
        FWHMeff            REAL,
        transparency       REAL,
        airmass            REAL,
        vSkyBright         REAL,
        filtSkyBright      REAL,
        rotSkyPos          REAL,
        rotTelPos          REAL,
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
        phaseAngle         REAL,
        rScatter           REAL,
        mieScatter         REAL,
        moonBright         REAL,
        darkBright         REAL,
        rawSeeing          REAL,
        wind               REAL,
        humidity           REAL,
        slewDist           REAL,
        slewTime           REAL,
        fiveSigmaDepth     REAL,
        ditheredRA         REAL,
        ditheredDec        REAL
        )""")

    def log_pointing(self, state, request):

        record = {}
        # TODO: might not want to use request_id here, but create a unique
        # non-null key
        record['obsHistID'] = request['request_id']
        record['sessionID'] = 0
        record['propID'] = request['target_program_id']
        record['fieldID'] = request['target_field_id']
        record['fieldRA'] = np.radians(request['target_ra'])
        record['fieldDec'] = np.radians(request['target_dec'])

        # TODO: double check if target_filter_id is letter or number
        record['filter'] = request['target_filter_id']
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

        pointing_seeing = seeing_at_pointing(state['current_zenith_seeing'].to(
            u.arcsec).value, altaz.alt.value)
        record['FWHMgeom'] = pointing_seeing
        record['FWHMeff'] = pointing_seeing
        # transparency

        # finRank
        record['airmass'] = altaz.secz.value
        # vSkyBright
        # record['filtSkyBrightness'] TODO
        record['rotSkyPos'] = 0.  # TODO: confirm
        record['rotTelPos'] = 0.
        record['lst'] = exposure_start.sidereal_time('apparent').to(
            u.hourangle).value
        record['altitude'] = altaz.alt.to(u.radian).value
        record['azimuth'] = altaz.az.to(u.radian).value

        sun = coord.get_sun(exposure_start)
        sun_altaz = skycoord_to_altaz(sun, exposure_start)
        moon = coord.get_moon(exposure_start, P48_loc)
        moon_altaz = skycoord_to_altaz(moon, exposure_start)
        record['dist2Moon'] = sc.separation(moon).to(u.radian).value
        record['solarElong'] = sc.separation(sun).to(u.deg).value
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
        # phaseAngle, rScatter, mieScatter, moonBright, darkBright
        # rawSeeing
        # wind
        # humidity
        if self.prev_obs is not None:
            sc_prev = coord.SkyCoord(self.prev_obs['fieldRA'] * u.radian,
                                     self.prev_obs['fieldDec'] * u.radian)
            record['slewDist'] = sc.separation(sc_prev).to(u.radian).value
            record['slewTime'] = (record['expDate'] -
                                  (self.prev_obs['expDate'] +
                                      self.prev_obs['visitTime']))
        # record['fiveSigmaDepth']
        record['ditheredRA'] = 0.
        record['ditheredDec'] = 0.

        # use placeholders to create the INSERT query
        columns = ', '.join(record.keys())
        placeholders = '{' + '}, {'.join(record.keys()) + '}'
        query = 'INSERT INTO Summary ({}) VALUES ({})'.format(
            columns, placeholders)
        query_filled = query.format(**record)
        self.conn.execute(query_filled)

        # save record for next obs
        self.prev_obs = record

        pass
