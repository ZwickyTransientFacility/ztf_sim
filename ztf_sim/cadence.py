
"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""

import numpy as np
from constants import *



def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True

def time_since_last_obs(request_row, current_state):
    """

    parameters for the cadence window are given in the cadence_pars
    dictionary of the row
    window_start : time delta (astropy quantity)
        time since last observation after which observations can begin
    window_stop : time delta (astropy quantity)
        time since last observation after which observations stop.
        Set to np.inf * u.days to have an open-ended window.
    prev_filter : {'same','any', 'other', [specific filterkey]}
        which filter to use to determine time of last observation
    """

    now = current_state['current_time']
    pars = row['cadence_pars']

    # which filter are we using to determine the time of last observation?
    # TODO: for now, require last observation to be from the same program
    if pars['prev_filter'] == 'any':
        last_obses = []
        for filter_id in FILTER_IDS:
            last_obses.append(row['last_observed_{}_{}'.format(
                row['program_id'],row['filter_id'])])
        last_obs = np.max(last_obses)
    else: 
        if pars['prev_filter'] == 'same':
            last_obs_filter = row['filter_id']
        elif pars['prev_filter'] == 'other':
            assert(len(FILTER_IDS) == 2)
            raise NotImplementedError
        elif pars['prev_filter'] in FILTER_IDS:
            last_obs_filter = pars['prev_filter']
   
        last_obs = row['last_observed_{}_{}'.format(
            row['program_id'],row['filter_id'])]

    window_start_ut = last_obs + pars['window_start']
    window_stop_ut = last_obs + pars['window_stop']

    return window_start_ut <= now <= window_stop_ut
