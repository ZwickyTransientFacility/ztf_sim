
import numpy as np
from sqlalchemy import create_engine
import sqlite3

class ObsLogger:

    def __init__(self, run_name, survey_start_time):
        self.run_name = run_name
        self.survey_start_time = survey_start_time
        self.prev_obs = None
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
        lst                REAL,
        alt                REAL,
        az                 REAL,
        dist2Moon          REAL,
        solarElong         REAL,
        moonRA             REAL,
        moonDec            REAL,
        moonAlt            REAL,
        moonAZ             REAL,
        moonPhase          REAL,
        sunAlt             REAL,
        sunAZ              REAL,
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
        record['fieldRA'] = request['target_ra']
        record['fieldDec'] = request['target_dec']
        # TODO: double check if target_filter_id is letter or number
        record['filter'] = request['target_filter_id']
        # times are recorded at start of exposure
        exposure_start = state['current_time'] - request['target_exposure_time']
        record['expDate'] = (exposure_start - self.survey_start_time).sec
        record['expMJD'] = exposure_start.mjd
        record['night'] = np.floor((exposure_start - self.survey_start_time).jd
                ).astype(np.int)
        record['visitTime'] = request['target_exposure_time'].value
        record['visitExpTime'] = request['target_exposure_time'].value
        # finRank
        record['FWHMgeom'] = state['current_seeing'].value
        record['FWHMeff'] = state['current_seeing'].value
        # transparency
        #record['airmass']

        # use placeholders to create the INSERT query
        columns = ', '.join(record.keys())
        placeholders = '{'+ '}, {'.join(record.keys()) + '}'
        query = 'INSERT INTO Summary ({}) VALUES ({})'.format(
                columns, placeholders)
        query_filled = query.format(**record)
        self.conn.execute(query_filled)

        pass
