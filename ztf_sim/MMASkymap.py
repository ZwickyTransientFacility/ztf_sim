"""MMA Skymaps."""

import os
from collections import defaultdict
from datetime import datetime
import logging
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time, TimeDelta
import astroplan
from .Fields import Fields

from .optimize import tsp_optimize, night_optimize
from .cadence import enough_gap_since_last_obs
from .constants import P48_loc, PROGRAM_IDS, FILTER_IDS, TIME_BLOCK_SIZE
from .constants import EXPOSURE_TIME, READOUT_TIME, FILTER_CHANGE_TIME, slew_time
from .constants import PROGRAM_BLOCK_SEQUENCE, LEN_BLOCK_SEQUENCE, MAX_AIRMASS
from .constants import BASE_DIR
from .utils import approx_hours_of_darkness
from .utils import skycoord_to_altaz, seeing_at_pointing
from .utils import altitude_to_airmass, airmass_to_altitude, RA_to_HA, HA_to_RA
from .utils import scalar_len, nightly_blocks, block_index, block_index_to_time
from .utils import block_use_fraction, maximum_altitude, compute_limiting_mag

class MMASkymap(object):

    def __init__(self, trigger_name, trigger_time, skymap_fields, fields=None):

        self.logger = logging.getLogger(__name__)

        self.trigger_name = trigger_name
        self.trigger_time = trigger_time
        self.skymap_fields = pd.DataFrame(skymap_fields)
        assert('field_id' in self.skymap_fields)
        assert('probability' in self.skymap_fields)

#        if fields is None:
#            self.fields = Fields()
#        else:
#            self.fields = fields

    def return_skymap(self):
        return self.skymap_fields

    def persist_skymap(self):
        pass
    
    def archive_persisted_skymap(self):
        pass



#def load_persisted_skymaps():
