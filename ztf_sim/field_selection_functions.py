"""Special-case field selection functions.

All functions should take 
    the current astropy.time.Time 
    an ObsLogger instance
    a dictionary of other observing programs field_ids and field_selection_functions
    a Fields object (for efficiency)
and return a list of field_ids."""

import logging
import numpy as np
import astropy.units as u
from astropy.time import Time
from .field_selection.srg import get_srg_fields
from .Fields import Fields
from .constants import PROGRAM_NAME_TO_ID, P48_Observer
from .utils import RA_to_HA

logger = logging.getLogger(__name__)

def msip_nss_selection_phaseii(time, obs_log, other_program_fields, fields):
    """Select MSIP NSS fields so we ensure lowdec coverage."""

    candidate_nss_field_ids = fields.select_field_ids(dec_range=[-32,90.],
                           grid_id=0,
                           # lowest rung fields are above airmass 2.5 for 2 hrs
                           observable_hours_range=[2.0, 24.])

    msip_internight_gap = 2*u.day
    msip_nobs_per_night = 2

    # unfortunate duplication of code from ObservingProgram.py
    # get the times they were last observed:
    # (note that fields *never* observed will not be included)
    # since this is function is for determining requests
    # at the start of the night, exclude observations taken tonight
    # this lets us restart the scheduler without breaking things
    last_observed_times = obs_log.select_last_observed_time_by_field(
        field_ids = candidate_nss_field_ids,
        filter_ids = [1,2],
        program_ids = [PROGRAM_NAME_TO_ID['MSIP']],
        subprogram_names = ['all_sky'],
        # arbitrary early date; start of night tonight
        mjd_range = [Time('2001-01-01').mjd,np.floor(time.mjd)])

    # we want an object observed at the end of the night N days ago
    # to be observed at the start of the night now.
    # Max night length is 12.2 hours
    cutoff_time = (time - (msip_internight_gap - 0.6 * u.day)).mjd

    # find fields last observed more recently than that
    wrecent = (last_observed_times['expMJD'] >= cutoff_time)
    recent_field_ids = last_observed_times.loc[wrecent].index.tolist()

    # reduce the list to only those not recently observed:
    nss_field_ids_due = [idi for idi in candidate_nss_field_ids if idi not in recent_field_ids] 
    nss_requests_allowed = other_program_fields[
            (PROGRAM_NAME_TO_ID['MSIP'],'all_sky')]['requests_allowed']
    n_fields_to_observe = int(np.round(nss_requests_allowed / msip_nobs_per_night))

    if len(nss_field_ids_due)*msip_nobs_per_night > nss_requests_allowed:
        # we have more fields we could observe than we are able to, select the
        # best ones
        logger.info(f'MSIP NSS: {nss_requests_allowed/msip_nobs_per_night} fields needed, {len(nss_field_ids_due)} available--removing fields.')
        available_nss_fields = fields.fields.loc[nss_field_ids_due]


        # find HA of the fields at midnight
        available_nss_fields['HA_midnight'] = RA_to_HA(
                available_nss_fields['ra'].values*u.degree, 
                P48_Observer.midnight(time, which='nearest')).wrap_at(180*u.deg)
        available_nss_fields['abs_HA_midnight'] = np.abs(available_nss_fields['HA_midnight'])
        nss_field_ids_to_observe = available_nss_fields.sort_values(by='abs_HA_midnight', ascending=True).iloc[:n_fields_to_observe].index.tolist()

        logger.info(f'MSIP NSS: requesting {nss_field_ids_to_observe}')

    else:
        # we have fewer fields available than expected--pad back on some recent
        # fields
        logger.info(f'MSIP NSS: {nss_requests_allowed/msip_nobs_per_night} fields needed, {len(nss_field_ids_due)} available--adding fields.')
        nss_field_ids_to_observe = nss_field_ids_due

        n_fields_needed = n_fields_to_observe - len(nss_field_ids_due)
        extra_fields = fields.fields.loc[recent_field_ids]

        # if we don't have enough extras it's okay; optimize.py will truncate
        # to the number of requests actually available
        if len(recent_field_ids) <= n_fields_needed:
            n_to_extend = len(recent_field_ids)
        else:
            n_to_extend = n_fields_needed

        # choose the ones with the lowest total observations?  Could make the
        # optimization tough since they will tend to only be up for short
        # periods
        nobs = obs_log.select_n_obs_by_field(field_ids=recent_field_ids,
                program_ids = [PROGRAM_NAME_TO_ID['MSIP']])
        extra_fields = extra_fields.join(nobs).fillna(0)

        field_ids_to_extend = extra_fields.sort_values(by='n_obs',ascending=True).iloc[:n_to_extend].index.tolist()
                
        logger.info(f'MSIP NSS: requesting {nss_field_ids_to_observe}')

        nss_field_ids_to_observe.extend(field_ids_to_extend)

        logger.info(f'MSIP NSS: Extending by {n_to_extend} fields: {field_ids_to_extend}')

    return nss_field_ids_to_observe

def msip_nss_selection(time, obs_log, other_program_fields, fields):
    """Select MSIP NSS fields so they don't overlap with other MSIP subprograms."""

    candidate_nss_fields = fields.select_field_ids(dec_range=[-32,90.],
                           grid_id=0,
                           observable_hours_range=[1.0, 24.])

    # now find the fields used by other MSIP subprograms
    other_MSIP_sp_fields = []
    for key, oops in other_program_fields.items():
        assert(len(key) == 2)
        if key[0] == PROGRAM_NAME_TO_ID['MSIP']:
            if key[1] != 'all_sky':
                if oops['field_ids'] is not None:
                    other_MSIP_sp_fields.extend(oops['field_ids'])
                else:
                    # we have to run the field selection function
                    # this is duplicative but easier than making a DAG
                    try:
                        selection_function = globals()[oops['field_selection_function']]
                        field_ids = selection_function(time, 
                            obs_log, other_program_fields, fields)
                        other_MSIP_sp_fields.extend(field_ids)
                    except Exception as e:
                        logger.exception(e)
                        logger.warning(f'During MSIP NSS field selection, error in generating nightly field list for {key}')

    nightly_nss_fields = np.setdiff1d(candidate_nss_fields, other_MSIP_sp_fields)

    logger.debug(f'MSIP NSS fields: {len(nightly_nss_fields)} are disjoint from other MSIP programs.')

    return nightly_nss_fields

def srg_selection(time, obs_log, other_program_fields, fields):
    """Select SRG fields."""

    # use the fields for multiple days so we can improve the sampling
    SRG_PAD_DAYS = 1
    dt = [0]
    for i in range(SRG_PAD_DAYS):
        dt.extend([i+1])
        dt.extend([-1*(i+1)])

    srg_fields = []
    for dti in dt:
        srg_fields.extend(get_srg_fields(time + dti * u.day, fields))
    
    srg_fields_unique = np.unique(srg_fields).tolist()

    logger.debug(f'SRG fields for {time.iso}+/-{SRG_PAD_DAYS} days: {srg_fields_unique}')

    return srg_fields_unique
