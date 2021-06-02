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

def msip_nss_selection_phaseii(time, obs_log, other_program_fields, fields,
        silent=False):
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
        if not silent:
            logger.info(f'MSIP NSS: {nss_requests_allowed/msip_nobs_per_night} fields needed, {len(nss_field_ids_due)} available--removing fields.')
        available_nss_fields = fields.fields.loc[nss_field_ids_due]
        available_nss_fields = available_nss_fields.join(last_observed_times)
        available_nss_fields['expMJD'] = available_nss_fields['expMJD'].fillna(Time('2001-01-01').mjd)
        nss_field_ids_to_observe = available_nss_fields.sort_values(by='expMJD',ascending=True).iloc[:n_fields_to_observe].index.tolist()

        if not silent:
            logger.info(f'MSIP NSS: requesting {nss_field_ids_to_observe}')

    else:
        # we have fewer fields available than expected--pad back on some recent
        # fields
        if not silent:
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
                
        if not silent:
            logger.info(f'MSIP NSS: requesting {nss_field_ids_to_observe}')

        nss_field_ids_to_observe.extend(field_ids_to_extend)

        if not silent:
            logger.info(f'MSIP NSS: Extending by {n_to_extend} fields: {field_ids_to_extend}')

    return nss_field_ids_to_observe

def partnership_HC_selection(time, obs_log, other_program_fields, fields):
    """Select partnership HC fields"""
    return phase_II_selection(time, obs_log, other_program_fields, fields,
                              subprogram='Partnership')

def Caltech_1DC_selection(time, obs_log, other_program_fields, fields):
    """Select Caltech 1DC fields"""
    return phase_II_selection(time, obs_log, other_program_fields, fields,
                              subprogram='Caltech')


def phase_II_selection(time, obs_log, other_program_fields, fields,
        subprogram='Partnership', silent=False):
    """Select partnership HC or Caltech 1DC fields"""

    assert (subprogram in ['Partnership', 'Caltech'])

    candidate_field_ids_prio1 = [473, 525, 624, 651, 671]
    candidate_field_ids_prio2 = [349, 350, 406, 414, 479, 480, 481, 528, 529, 530, 531, 532, 533, 534, 544,
545, 546, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 596, 597, 600, 601,
623, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 644, 645, 646, 647,
668, 669, 670, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 708, 709, 710,
711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 748, 749,
750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 784, 785, 786,
787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 815, 816, 817, 818, 819, 820,
821, 822, 823, 824, 825, 841, 843, 844]
    candidate_field_ids_prio3 = [345, 346, 347, 348, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401,
402, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 441, 442,
443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 464, 465, 466, 467, 468, 469,
470, 471, 472, 474, 475, 476, 477, 478, 493, 494, 496, 497, 498, 499, 500, 501,
502, 503, 516, 517, 518, 519, 520, 521, 522, 523, 524, 526, 527, 547, 548, 549,
550, 551, 552, 553, 554, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 599,
602, 604, 615, 616, 617, 618, 619, 620, 621, 622, 662, 663, 664, 665, 666, 667,
706, 707]
    candidate_field_ids_prio4 = [296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 351, 352, 353, 354,
355, 403, 405, 413, 482, 483, 484, 535, 543, 586, 594, 595, 682, 746, 747, 763,
782, 783, 814, 839, 840, 842, 845, 846, 847, 848, 858, 859, 860, 861, 862, 863,
864]

    candidate_field_ids_byprio = [candidate_field_ids_prio1, candidate_field_ids_prio2, candidate_field_ids_prio3, candidate_field_ids_prio4]

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

    # define the parameters that vary
    if subprogram == 'Partnership':
        label = 'Partnership HC'
        min_hours_visible = 2.5
        nobs_per_night = 4
        requests_allowed = other_program_fields[
            (PROGRAM_NAME_TO_ID['collaboration'],'high_cadence')]['requests_allowed']
        offset = 0
        exclude_fields = []


    elif subprogram == 'Caltech':
        label = 'Caltech 1DC'
        min_hours_visible = 1.5
        nobs_per_night = 2
        requests_allowed = other_program_fields[
            (PROGRAM_NAME_TO_ID['Caltech'],'Caltech_1DC')]['requests_allowed']
        # to prevent fields from thrashing back and forth between programs
        # we place a buffer in the priority list
        offset = 15
        exclude_fields = phase_II_selection(time, obs_log, other_program_fields, fields, subprogram='Partnership', silent=False)


    else:
        raise ValueError(f'Unknown subprogram {subprogram}')


    # do a visibility check (also handles moon exclusion)
    visible_field_ids = fields.select_field_ids(dec_range=[-32,90.],
                           grid_id=0,
                           observable_hours_range=[min_hours_visible, 24.])

    visible_field_ids = np.intersect1d(visible_field_ids, 
            candidate_field_ids).tolist()

    n_fields_to_observe = int(np.round(requests_allowed / nobs_per_night))
    n_fields_remaining = n_fields_to_observe

    if not silent:
        logger.info(f'{label}: need {n_fields_to_observe} fields')

    field_ids_to_observe = []

    for i, prio_field_ids in enumerate(candidate_field_ids_byprio):

        assert (n_fields_remaining >= 0)
        if n_fields_remaining == 0:
            break

        visible_prio_field_ids = np.intersect1d(visible_field_ids,
            prio_field_ids).tolist()

        # for Caltech, exclude partnership fields being observed
        visible_prio_field_ids = [f for f in visible_prio_field_ids if f not in exclude_fields]

        if len(visible_prio_field_ids) == 0:
            if not silent:
                logger.info(f'{label}: No fields remaining at priority {i+1}')
            continue

        if len(visible_prio_field_ids) <= offset:
            #always observe prio any 1 fields if available, without offset
            # they weren't included in Partnership HC due to min_hours_visible
            if i==0: # prio 1
                field_ids_to_observe.extend(visible_prio_field_ids)
                n_fields_remaining -= len(visible_prio_field_ids)
                if not silent:
                    logger.info(f'{label}: requesting {len(visible_prio_field_ids)} priority {i+1} fields: {visible_prio_field_ids}')

            else:
                offset -= len(visible_prio_field_ids)
                if not silent:
                    logger.info(f'{label}: Offsetting {len(visible_prio_field_ids)} at priority {i+1}')

            continue

        if (len(visible_prio_field_ids)+offset) <= n_fields_remaining:
            # include all fields in this priority
            field_ids_to_observe.extend(visible_prio_field_ids[offset:])
            n_fields_remaining -= len(visible_prio_field_ids[offset:])
            if not silent and offset != 0:
                logger.info(f'{label}: Offsetting {offset} at priority {i+1}')
            offset = 0
            if not silent:
                logger.info(f'{label}: requesting {len(visible_prio_field_ids[offset:])} priority {i+1} fields: {visible_prio_field_ids}')
        else:
            # prioritize by cutting on HA at midnight
            
            if not silent and offset != 0:
                logger.info(f'{label}: Offsetting {offset} at priority {i+1}')

            candidate_fields = fields.fields.loc[visible_prio_field_ids].copy()

            candidate_fields['HA_midnight'] = RA_to_HA(
                        candidate_fields['ra'].values*u.degree, 
                    P48_Observer.midnight(time, which='nearest')).wrap_at(180*u.deg)
            candidate_fields['abs_HA_midnight'] = np.abs(candidate_fields['HA_midnight'])
            footprint_field_ids = candidate_fields.sort_values(
                    by='abs_HA_midnight', ascending=True).iloc[offset:offset+n_fields_remaining].index.tolist()
            offset=0
            n_fields_remaining -= len(footprint_field_ids)
            field_ids_to_observe.extend(footprint_field_ids)
            if not silent:
                logger.info(f'{label}: requesting {len(footprint_field_ids)} priority {i+1} fields of {len(candidate_fields)} total: {footprint_field_ids}')

    if not silent:
        logger.info(f'{label}: requesting {field_ids_to_observe}')

    return field_ids_to_observe

def srg_selection(time, obs_log, other_program_fields, fields):
    """Select SRG fields."""

    # use the fields for multiple days so we can improve the sampling
    SRG_PAD_DAYS = 2
    dt = [0]
    for i in range(SRG_PAD_DAYS):
        dt.extend([i+1])
        dt.extend([-1*(i+1)])

    srg_fields = []
    for dti in dt:
        srg_fields.extend(get_srg_fields(time + dti * u.day, fields))
    
    srg_fields_unique = np.unique(srg_fields).tolist()

    # exclude fields being observed by MSIP
    nss_selection_function = globals()[other_program_fields[(PROGRAM_NAME_TO_ID['MSIP'],'all_sky')]['field_selection_function']]
    nss_field_ids = nss_selection_function(time, 
        obs_log, other_program_fields, fields, silent=True)

    srg_fields = np.setdiff1d(srg_fields_unique, nss_field_ids)

    logger.debug(f'SRG fields for {time.iso}+/-{SRG_PAD_DAYS} days: {srg_fields}')

    return srg_fields


def aam_caltech_june21(time, obs_log, other_program_fields, fields):
    """Select Ashish fields, June 2021."""

    aam_fields = []

    # determine fields being observed by MSIP
    nss_selection_function = globals()[other_program_fields[(PROGRAM_NAME_TO_ID['MSIP'],'all_sky')]['field_selection_function']]
    nss_field_ids = nss_selection_function(time, 
        obs_log, other_program_fields, fields, silent=True)

    if (time >= Time('2021-06-05')) and (time < Time('2021-06-11')):
        # Daily from Jun 4 through Jun 9 PT inclusive (== Jun 5 through Jun 10 UT - 6 nights) 3 ZTF fields viz. 1433, 1476, and 1525 once each in g and r (g and r need not be consecutive, but close in time will be useful).
        # set 1 fields: observe every night
        set1 = [1433, 1476, 1525]
        aam_fields.extend(set1)

        # Fields 530, 482, 388 from Jun 4 through Jun 9 PT (== Jun 5 through Jun 10 UT, 6 nights), one image each field in g and r on the days when these fields are not being covered by the MSIP survey (the alternate days). Again, g and r need not be consecutive in time.

        # set 2 fields: observe when not being observed by MSIP
        set2 = [530, 482, 388]
        set2_use = np.setdiff1d(set2, nss_field_ids)

        aam_fields.extend(set2_use)

    if (time >= Time('2021-06-03')) and (time < Time('2021-06-05')):
        #(3) Fields 1433, 1476, and 1525 once in g and r on Jun 2 PT (Jun 3 UT)
        #If these can not be taken on Jun 2 PT, we request that the observations are attempted on Jun 3 PT (Jun 4 UT).
        set3 = [1433, 1476, 1525]

        observed_ids = obs_log.select_last_observed_time_by_field(
            field_ids = set3,
            filter_ids = [1,2],
            program_ids = [PROGRAM_NAME_TO_ID['Caltech']],
            subprogram_names = ['AAM_June21'],
            # arbitrary early date; start of night tonight
            mjd_range = [Time('2001-01-01').mjd,np.floor(time.mjd)]).index.tolist()

        set3_use = np.setdiff1d(set3, observed_ids)

        aam_fields.extend(set3_use)

    if (time >= Time('2021-06-11')) and (time < Time('2021-07-16')):
        #(4) Fields 1433, 1476, and 1525 in g and r once per week for five weeks after Jun 10 (ideally on days when Fields 530, 482, 388 are not done). 
        set4 = [1433, 1476, 1525]
        set4_use = set4

        # check if it's been observed in the last week under this program

        # unfortunate duplication of code from ObservingProgram.py
        # get the times they were last observed:
        # (note that fields *never* observed will not be included)
        # since this is function is for determining requests
        # at the start of the night, exclude observations taken tonight
        # this lets us restart the scheduler without breaking things
        last_observed_times = obs_log.select_last_observed_time_by_field(
            field_ids = set4_use,
            filter_ids = [1,2],
            program_ids = [PROGRAM_NAME_TO_ID['Caltech']],
            subprogram_names = ['AAM_June21'],
            # arbitrary early date; start of night tonight
            mjd_range = [Time('2001-01-01').mjd,np.floor(time.mjd)])

        # we want an object observed at the end of the night N days ago
        # to be observed at the start of the night now.
        # Max night length is 12.2 hours
        intranight_gap = 7. * u.day
        cutoff_time = (time - (intranight_gap - 0.6 * u.day)).mjd

        # find fields last observed more recently than that
        wrecent = (last_observed_times['expMJD'] >= cutoff_time)
        recent_field_ids = last_observed_times.loc[wrecent].index.tolist()

        # reduce the list to only those not recently observed:
        field_ids_due = [idi for idi in set4_use if idi not in recent_field_ids] 
        aam_fields.extend(field_ids_due)

        logger.info(f'Caltech AAM: requesting {aam_fields}')

    return aam_fields




##################### previous/unused -- for reference

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
