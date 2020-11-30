"""Special-case field selection functions.

All functions should take 
    the current astropy.time.Time 
    an ObsLogger instance
    a dictionary of other observing programs field_ids and field_selection_functions
    a Fields object (for efficiency)
and return a list of field_ids.

Note that any cadence cuts should occur within the functions defined here--the cadence cuts in ObservingProgram.py are overridden. 
"""

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
    candidate_nss_fields = fields.fields.loc[candidate_nss_field_ids].copy()


    msip_cadence = 2
    msip_internight_gap = msip_cadence*u.day
    msip_nobs_per_night = 2

    nss_requests_allowed = other_program_fields[
            (PROGRAM_NAME_TO_ID['MSIP'],'all_sky')]['requests_allowed']
    n_fields_to_observe = int(np.round(nss_requests_allowed / msip_nobs_per_night))

    # define a footprint for the 2-night cadence by cutting on HA at midnight
    candidate_nss_fields['HA_midnight'] = RA_to_HA(
                candidate_nss_fields['ra'].values*u.degree, 
            P48_Observer.midnight(time, which='nearest')).wrap_at(180*u.deg)
    candidate_nss_fields['abs_HA_midnight'] = np.abs(candidate_nss_fields['HA_midnight'])
    nss_footprint_field_ids = candidate_nss_fields.sort_values(by='abs_HA_midnight', ascending=True).iloc[:(msip_cadence*n_fields_to_observe-1)].index.tolist()
    # save these if we need filler
    candidate_field_ids_outside_footprint = candidate_nss_fields.sort_values(by='abs_HA_midnight', ascending=True).iloc[(msip_cadence*n_fields_to_observe-1):].index.tolist()

    # unfortunate duplication of code from ObservingProgram.py
    # get the times they were last observed:
    # (note that fields *never* observed will not be included)
    # since this is function is for determining requests
    # at the start of the night, exclude observations taken tonight
    # this lets us restart the scheduler without breaking things
    last_observed_times = obs_log.select_last_observed_time_by_field(
        field_ids = nss_footprint_field_ids,
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
    nss_field_ids_due = [idi for idi in nss_footprint_field_ids if idi not in recent_field_ids] 

    if len(nss_field_ids_due)*msip_nobs_per_night > nss_requests_allowed:
        # we have more fields we could observe than we are able to
        # select the fields with the least recent observations
        logger.info(f'MSIP NSS: {nss_requests_allowed/msip_nobs_per_night} fields needed, {len(nss_field_ids_due)} available--removing fields.')
        available_nss_fields = fields.fields.loc[nss_field_ids_due]
        available_nss_fields = available_nss_fields.join(last_observed_times)
        available_nss_fields['expMJD'] = available_nss_fields['expMJD'].fillna(Time('2001-01-01').mjd)
        nss_field_ids_to_observe = available_nss_fields.sort_values(by='expMJD',ascending=True).iloc[:n_fields_to_observe].index.tolist()

        logger.info(f'MSIP NSS: requesting {nss_field_ids_to_observe}')

    else:
        # we have fewer fields available than expected--pad back on some recent
        # fields
        logger.info(f'MSIP NSS: {nss_requests_allowed/msip_nobs_per_night} fields needed, {len(nss_field_ids_due)} available--adding fields.')
        nss_field_ids_to_observe = nss_field_ids_due

        n_fields_needed = n_fields_to_observe - len(nss_field_ids_due)

        # for now just pad on the extra fields outside the footprint
        # TODO: think more carefully about how best to handle this
        recent_field_ids.extend(candidate_field_ids_outside_footprint)

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

def partnership_HC_selection(time, obs_log, other_program_fields, fields):
    """Select partnership HC fields"""

    candidate_field_ids = [296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 345, 346, 347, 348, 349,
350, 351, 352, 353, 354, 355, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400,
401, 402, 403, 405, 406, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423,
424, 425, 426, 427, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452,
464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479,
480, 481, 482, 483, 484, 493, 494, 496, 497, 498, 499, 500, 501, 502, 503, 516,
517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532,
533, 534, 535, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 566,
567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582,
583, 584, 585, 586, 594, 595, 596, 597, 599, 600, 601, 602, 604, 615, 616, 617,
618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633,
634, 635, 644, 645, 646, 647, 651, 662, 663, 664, 665, 666, 667, 668, 669, 670,
671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 706, 707, 708, 709,
710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 746,
747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762,
763, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796,
814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 839, 840, 841, 842,
843, 844, 845, 846, 847, 848, 858, 859, 860, 861, 862, 863, 864]
    candidate_fields = fields.fields.loc[candidate_nss_field_ids].copy()


    HC_cadence = 1
    HC_internight_gap = HC_cadence*u.day
    HC_nobs_per_night = 4

    HC_requests_allowed = other_program_fields[
            (PROGRAM_NAME_TO_ID['collaboration'],'high_cadence')]['requests_allowed']
    n_fields_to_observe = int(np.round(HC_requests_allowed / HC_nobs_per_night))

    # define a footprint for the 2-night cadence by cutting on HA at midnight
    candidate_fields['HA_midnight'] = RA_to_HA(
                candidate_fields['ra'].values*u.degree, 
            P48_Observer.midnight(time, which='nearest')).wrap_at(180*u.deg)
    candidate_fields['abs_HA_midnight'] = np.abs(candidate_fields['HA_midnight'])
    HC_footprint_field_ids = candidate_fields.sort_values(by='abs_HA_midnight', ascending=True).iloc[:(HC_cadence*n_fields_to_observe-1)].index.tolist()

    logger.info(f'Partnership HC: requesting {HC_footprint_field_ids}')

    return HC_footprint_field_ids

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
