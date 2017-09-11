
"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""
from __future__ import absolute_import

import numpy as np
from .constants import FILTER_IDS


def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True


def time_since_obs(request_row, current_state):
    """

    parameters for the cadence window are given in the cadence_pars
    dictionary of the row:
        'ref_obs': ['last_observed', 'first_obs_tonight']
            which previous observation to use as reference time
        'prev_filter': ['same','other', 'any']
            which filter the previous observation should be in
        'window_start', 'window_stop'
            relative time bounds for the window after the ref_obs time
            e.g, window_start = 1.0*u.hour starts at ref_obs + 1 hour
        'filter_id'
            which filter to observe in
        'program_id'
            program requesting the observation
    One pair of :
        window_start : time delta (astropy quantity)
            time since ref observation after which observations can begin
        window_stop : time delta (astropy quantity)
            time since ref observation after which observations stop.
    or:
        window_center : time (astropy quantity)
            middle of acceptable window since ref observation
        window_half_width : time delta (astropy quantity)
            half width of acceptable window.

    prev_filter : {'same','any', 'other', [specific filterkey]}
        which filter to use to determine time of last observation
    """

    now = current_state['current_time'].mjd
    pars = request_row['cadence_pars']

    ref_obs = get_ref_obs_time(request_row, current_state)

    if np.isnan(ref_obs):
        # for first_obs_tonight scheduling, have to wait until first obs is
        # taken
        return False

    if 'window_start' in pars and 'window_stop' in pars:
        window_start_ut = ref_obs + pars['window_start']
        window_stop_ut = ref_obs + pars['window_stop']
    elif 'window_center' in pars and 'window_half_width' in pars:
        window_start_ut = ref_obs + \
            pars['window_center'] - pars['window_half_width']
        window_stop_ut = ref_obs + \
            pars['window_center'] + pars['window_half_width']

    return window_start_ut <= now <= window_stop_ut

def get_ref_obs_time(request_row, current_state):
    """Determine which past observation should be used to make the cadence window reference time."""

    # which filter are we using to determine the time of last observation?
    # TODO: for now, require last observation to be from the same program
    pars = request_row['cadence_pars']

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
