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
from .field_selection.srg import get_srg_fields
from .Fields import Fields
from .constants import PROGRAM_NAME_TO_ID

logger = logging.getLogger(__name__)

def msip_nss_selection(time, obs_log, other_program_fields, fields):
    """Select MSIP NSS fields so they don't overlap with other MSIP subprograms."""

    candidate_nss_fields = fields.select_field_ids(dec_range=[-31.,90.],
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
