"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""

import pandas as pd
import numpy as np
import astropy.units as u
from .constants import FILTER_IDS, TIME_BLOCK_SIZE


def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True


def enough_gap_since_last_obs(df, current_state, obs_log):
    """
    Determine if a sufficient time has passed since the last observation
    in this subprogram (in any filter):
    """

    now = current_state['current_time'].mjd

    # don't mess up with the upstream data structure
    df = df.copy()

    grp = df.groupby(['program_id','subprogram_name'])

    df['ref_obs_mjd'] = np.nan
    for grpi, dfi in grp:
        ref_obs = obs_log.select_last_observed_time_by_field(
                field_ids = set(dfi['field_id'].tolist()), 
                program_ids = [grpi[0]],
                subprogram_names = [grpi[1]])
        if len(ref_obs) > 0:
            tmp = pd.merge(df, ref_obs, left_on='field_id', right_index=True,
                    how='inner')
            df.loc[tmp.index, 'ref_obs_mjd'] = tmp.expMJD.values


    # give a fake value for fields unobserved
    df.loc[df['ref_obs_mjd'].isnull(), 'ref_obs_mjd'] = 58119.0

    # calculate dt
    df['dt'] = now - df['ref_obs_mjd']

    return df['dt'] >=  (df['intranight_gap_min']*(1*u.minute).to(u.day).value)
