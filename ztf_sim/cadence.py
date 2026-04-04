"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""

import pandas as pd
import numpy as np
import astropy.units as u
from .constants import FILTER_IDS, TIME_BLOCK_SIZE


def no_cadence(*args):
    """Return True unconditionally, indicating no cadence constraint.

    Drop-in replacement for `enough_gap_since_last_obs` for programs that
    have no intranight gap requirement.

    Parameters
    ----------
    *args
        Accepted but ignored, for signature compatibility.

    Returns
    -------
    bool
        Always ``True``.
    """
    return True


def enough_gap_since_last_obs(df, current_state, obs_log):
    """Check whether sufficient intranight time has elapsed since each field's last observation.

    Groups candidates by ``(program_id, subprogram_name)`` and queries the
    observation log for the most recent exposure of each field within the
    current night. Fields that have never been observed tonight are always
    treated as eligible.

    Parameters
    ----------
    df : pandas.DataFrame
        Candidate observations. Required columns: ``program_id`` (int),
        ``subprogram_name`` (str), ``field_id`` (int),
        ``intranight_gap_min`` (float, minutes).
    current_state : dict
        Telescope state dict as returned by
        ``TelescopeStateMachine.current_state_dict()``. Must contain key
        ``'current_time'`` (astropy.time.Time).
    obs_log : ObsLogger
        Observation history used to look up the last observed time per field.

    Returns
    -------
    pandas.Series of bool
        Boolean Series with the same index as *df*. ``True`` where the
        elapsed time since the last observation meets or exceeds
        ``intranight_gap_min``.
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
