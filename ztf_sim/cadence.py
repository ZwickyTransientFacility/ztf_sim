
"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""
from __future__ import absolute_import

import numpy as np
import astropy.units as u
from .constants import FILTER_IDS, TIME_BLOCK_SIZE


def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True


def enough_gap_since_last_obs(request_row, current_state):
    """
    Determine if a sufficient time has passed since the last observation
    in this program (in any filter):
        'filter_id'
            which filter to observe in
        'program_id'
            program requesting the observation
    """

    now = current_state['current_time'].mjd
    pars = {'prev_filter':'any','ref_obs':'last_observed',
            'min_gap':TIME_BLOCK_SIZE}

    ref_obs = get_ref_obs_time(request_row, current_state, pars)

    if np.isnan(ref_obs):
        # for first_obs_tonight scheduling (currently unimplemented), 
        # have to wait until first obs is taken
        return False

    return now >= (ref_obs + pars['min_gap'].to(u.day).value)

def get_ref_obs_time(request_row, current_state, pars):
    """Determine which past observation should be used to make the cadence window reference time."""

    if pars['prev_filter'] == 'any':
        ref_obses = []
        for filter_id in FILTER_IDS:
            ref_obses.append(request_row['{}_{}_{}'.format(
                pars['ref_obs'],
                request_row['program_id'], filter_id)])
        ref_obses = np.array(ref_obses)
        if np.sum(np.isfinite(ref_obses)):
            ref_obs = np.max(ref_obses[np.isfinite(ref_obses)])
        else:
            # no previous observations
            ref_obs = np.nan
    else:
        if pars['prev_filter'] == 'same':
            ref_obs_filter = request_row['filter_id']
        elif pars['prev_filter'] == 'other':
            assert(len(FILTER_IDS) == 2)
            raise NotImplementedError
        elif pars['prev_filter'] in FILTER_IDS:
            ref_obs_filter = pars['prev_filter']

        ref_obs = request_row['{}_{}_{}'.format(
            pars['ref_obs'],
            request_row['program_id'], ref_obs_filter)]

    return ref_obs

def absolute_time_window(request_row, current_state):
    """stub for assigning fields specific UTC slots"""
    raise NotImplementedError
