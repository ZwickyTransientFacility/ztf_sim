
"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""
from __future__ import absolute_import

import numpy as np
import astropy.units as u
from .constants import FILTER_IDS, TIME_BLOCK_SIZE


def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True


def enough_gap_since_last_obs(request_row, current_state, obs_log):
    """
    Determine if a sufficient time has passed since the last observation
    in this subprogram (in any filter):
    """

    now = current_state['current_time'].mjd
    min_gap = TIME_BLOCK_SIZE

    ref_obs = obs_log.select_last_observed_time_by_field(
            field_ids = [request_row['field_id']], 
            program_ids = [request_row['program_id']],
            subprogram_names = [request_row['subprogram_name']])

    if len(ref_obs) == 0:
        # field has never been observed in this subprogram, 
        # it's okay by definition
        return True

    assert( (len(ref_obs) == 1) )

    ref_obs = ref_obs.expMJD.values[0]

    return now >= (ref_obs + min_gap.to(u.day).value)
