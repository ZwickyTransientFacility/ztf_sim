import numpy as np
import glob, re
import pandas as pd
import requests
import os
import time
from scipy.stats import mode

from gurobipy import *


from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.data import get_readable_fileobj
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun

from astroplan import FixedTarget, Observer, Constraint, is_observable, min_best_rescale,AirmassConstraint

from ..Fields import Fields
from ..constants import READOUT_TIME
from ..configuration import Configuration
from ..optimize import check_limits_and_solve_TSP
from ..QueueManager import ListQueueManager

logger = logging.getLogger(__name__)



def make_ep_blocks(time_now, time_allowed, time_limit=300*u.second,
                   other_timed_queues_tonight = []):

    ##############################

    # getting ZTF fields

    f = Fields()

    cuts = f.select_fields(grid_id=0, dec_range=[-30, 90])
    ztf_fields = f.fields[cuts]
    # index of ztf_fields are the field_ids matching coord_ztf_centers
    coord_ztf_centers = f._field_coords(cuts=cuts)

    # Download the EP schedule
    EP_schedule = 'https://ep.bao.ac.cn/ep/observation_plan/download_obsplan_fov'
    df_ep_fov = None
    for i in range(5):
        try:
            logging.info(f"Attempting EP download, attempt {i}")
            df_ep_fov = pd.read_json(EP_schedule)
        except Exception as e:
            logging.exception(e)
            time.sleep(30)
        else:
            break

    if df_ep_fov is None:
        raise ValueError("EP pointings failed to download!")

    logging.info(f"Downloaded {len(df_ep_fov)} EP pointings.")
    #df_ep_fov.to_json(f'../../sims/EP_{time_now.isot}.json')
    # provide a unique index for the pointing
    df_ep_fov['ep_pointing_id'] = df_ep_fov.index
    df_ep_fov['start_mjd'] = [Time(d).mjd for d in df_ep_fov['start_date_UTC']]
    df_ep_fov['end_mjd'] = [Time(d).mjd for d in df_ep_fov['end_date_UTC']]
    df_ep_fov['tobs_sec'] = (df_ep_fov['end_mjd'] - df_ep_fov['start_mjd']) * 24* 3600.

    palomar = Observer.at_site('palomar')
    dayfrac = time_now.mjd % np.floor(time_now.mjd)
    # Palomar sunset is after zero UTC 
    time_zero_utc = Time(np.floor(time_now.mjd), format='mjd')
    if dayfrac > 0.5: 
        # we're in the daytime, try for tonight 
        time_zero_utc +=1

    sunset_18 = palomar.sun_set_time(time_zero_utc, which='next', horizon=-18*u.deg)
    sunrise_18 = palomar.sun_rise_time(time_zero_utc, which='next', horizon=-18*u.deg)
    airmass = 2

    cond_sunset = df_ep_fov.start_date_UTC > sunset_18
    cond_sunrise = df_ep_fov.end_date_UTC < sunrise_18
    cond_night = cond_sunrise &  cond_sunset
    logging.info(f"EP: {np.sum(cond_night)} nighttime pointings.")
    logging.info(f"EP: {df_ep_fov.loc[cond_night, ['start_mjd', 'end_mjd']]}")

    # check for interference with other timed queues

    if len(other_timed_queues_tonight) > 0:
        # make a set of Trues
        cond_no_intersection = df_ep_fov['ep_pointing_id'] == df_ep_fov['ep_pointing_id']
        for oq in other_timed_queues_tonight:
            # Neutrino TOOs tend to be short but have long validity windows
            # I could compute the exposure time in the queue and do a cut, 
            # but I think I'll just exclude all ToOs since they will preempt EP anyway
            if oq.queue_name.startswith('ToO'):
                continue

            # https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
            # timed queue doesn't intersect EP pointing if it ends before the EP 
            # EP pointing starts or starts after the EP pointing ends
            logging.info(f"{oq.queue_name}, {oq.validity_window[0].mjd} - {oq.validity_window[1].mjd}")
            cond_no_intersection &= ((oq.validity_window[1].mjd < df_ep_fov['start_mjd']) | 
                                     (df_ep_fov['end_mjd'] < oq.validity_window[0].mjd))

        logging.info(f"{np.sum(~cond_no_intersection & cond_night)} potential nighttime EP pointings intersect with timed queues.")
        cond_good = cond_night & cond_no_intersection

    else:
        cond_good = cond_night

    df_ep_fov = df_ep_fov[cond_good]
    ep_pointing_idxs = df_ep_fov.index.tolist()

    logging.info(f"Looking for overlap with {len(ep_pointing_idxs)} EP pointings.")

    if len(ep_pointing_idxs) == 0:
        raise ValueError("No available EP pointings!")

    # for each EP pointing:
    #     get which ZTF fields are observable
    #     save field, time

    matched_dfs = []

    for wxt_pointing in ep_pointing_idxs:
        
        ep_fov_i = np.array(df_ep_fov.loc[wxt_pointing,'WXT_fov'])

        coord_wxt = SkyCoord(ep_fov_i, unit='deg')

        # find the closest ZTF field to each EP fov point
        idxz, d2d, _ = coord_wxt.match_to_catalog_sky(coord_ztf_centers)
        wmatch = d2d < 3*u.degree
        # store the unique matches (index into coord_ztf_centers)
        uniq_match_idxz = np.unique(idxz[wmatch])

        # if no overlap continue
        if np.sum(wmatch) == 0:
           continue 
            

        # check if the ZTF fields are observable
        
        t_start = Time(df_ep_fov.start_date_UTC[wxt_pointing])
        t_end = Time(df_ep_fov.end_date_UTC[wxt_pointing])
        

        observable = is_observable(AirmassConstraint(airmass),palomar,
                                   coord_ztf_centers[uniq_match_idxz], t_start)

        if np.sum(observable) == 0:
            continue

        # observable field ids
        zfi = ztf_fields.index[uniq_match_idxz].values[observable].tolist()

        # count the number of EP FOV points for each field
        ep_fov_coverage = []
        for i, zf in enumerate(uniq_match_idxz):
            if observable[i]:
                separations = coord_ztf_centers[zf].separation(coord_wxt)
                ep_fov_coverage.append(np.sum(separations < 3.75*u.degree))
        
        assert (len(zfi) == len(ep_fov_coverage))




        # store the results in a tidy data frame

        dict_ztf_wxt_temp = {
                'field_id': zfi,
                'ep_pointing_id': np.zeros(len(zfi), dtype=int) + wxt_pointing,
                'ep_fov_coverage': ep_fov_coverage,
                'ra': ztf_fields.loc[zfi, 'ra'].tolist(),
                'dec': ztf_fields.loc[zfi, 'dec'].tolist()
                }

        matched_dfs.append(pd.DataFrame(dict_ztf_wxt_temp))

        
    dfm = pd.concat(matched_dfs, ignore_index=True)

    unique_ztf_fields = dfm['field_id'].unique()


    # Create an empty model
    m = Model('ep')

    # set the number of threads Gurobi uses
    if 'GRB_USE_NTHREADS' in os.environ:
        m.Params.Threads = int(os.environ['GRB_USE_NTHREADS'])

    # decision variable: yes or no for each field f in pointing i
    # note that I'm being a little sloppy here: the variable name here
    # is the index of dfm, which is not actually a unique field id
    yfi_dict = m.addVars(dfm.index,name='Yfi',vtype=GRB.BINARY)
    yfi_series = pd.Series(yfi_dict,name='Yfi')
    dfm = dfm.join(yfi_series)

    # putting the Yfi vars in the dfm series is not always convenient;
    # provide a lookup table back to the dfm index
    fi_to_idx = {}
    idx_to_fi = dfm[['field_id','ep_pointing_id']].to_dict(orient='index')
    for k,v in idx_to_fi.items():
        fi_to_idx[(v['field_id'],v['ep_pointing_id'])] = k

    # yes or no if each pointing i gets ZTF observations
    yi_dict = m.addVars(ep_pointing_idxs,name='Yi',vtype=GRB.BINARY)
    yi_series = pd.Series(yi_dict,name='Yi', index=ep_pointing_idxs)
    df_ep_fov = df_ep_fov.join(yi_series)


    # create resultant variables: Yi = 1 if ep pointing i has 
    # at least one observation
    for i in ep_pointing_idxs:
        fi_set = dfm.loc[dfm['ep_pointing_id'] == i, 'field_id'].values.tolist()
        for f in fi_set:
            m.addGenConstrOr(yi_dict[i], dfm.loc[(dfm['ep_pointing_id'] == i) & 
                                                 (dfm['field_id'] == f), 'Yfi'],
                "orconstr_{}".format(i))

    # create resultant variables: 1 if field if is observed in multiple pointings
    yfik = m.addVars(unique_ztf_fields,ep_pointing_idxs[:-1],ep_pointing_idxs[1:],
            vtype=GRB.BINARY)
    constraint_dict = {}
    for f in unique_ztf_fields:
        for i in ep_pointing_idxs[:-1]:
            for k in ep_pointing_idxs[1:]:
                if i < k:
                    # avoid duplicate entries
                    continue
                try:
                    m.addGenConstrAnd(yfik[f,i,k],
                        [yfi_dict[fi_to_idx[(f,i)]],
                         yfi_dict[fi_to_idx[(f,k)]]],
                        "jointobs_and_{}_{}_{}".format(f,i,k))
                except KeyError:
                    # not all fields are in all ep pointings
                    continue

    # total exposure time constraint
    constr_total_exp = m.addConstr(
        (np.sum(dfm['Yfi'] *
            (30 + READOUT_TIME.to(u.second).value))
            <= (time_allowed.to(u.second).value)) , "constr_total_exp")

    # don't overfill each EP pointing
    tobs_pointing = df_ep_fov.loc[:,'tobs_sec'].to_dict()
    constr_nperslot = m.addConstrs(
        (np.sum(dfm.loc[dfm['ep_pointing_id'] == i, 'Yfi'] *
            (30 + READOUT_TIME.to(u.second).value))
            <= tobs_pointing[i]
            for i in ep_pointing_idxs), "constr_nperslot")

    # don't make too many blocks
    # pointings seem to be about 1500-4k seconds
    # 5-10 pointings would mean 20% slack would give 1-2 extra blocks
    block_slack = 1.2
    constr_n_ep = m.addConstr(
        (np.sum(df_ep_fov['Yi']*df_ep_fov['tobs_sec']) <=
            (block_slack * time_allowed.to(u.second).value)) , "constr_n_exp")


    
    m.update()



    # ep_fov_coverage ranges from 0-56 or so
    repeat_bonus = 50
    m.setObjective(
            np.sum(dfm['Yfi'] * dfm['ep_fov_coverage']) 
            + repeat_bonus * yfik.sum(),
            GRB.MAXIMIZE)

    m.Params.TimeLimit = time_limit.to(u.second).value

    m.update()

    m.optimize()

    if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
        logger.warning(f"Gurobi failed to optimize! Code {m.Status}")

    # now get the decision variables.  Use > a constant to avoid
    # numerical precision issues
    dfm['Yfi_val'] = dfm['Yfi'].apply(lambda x: x.getAttr('x') > 0.1)

    df_schedule = dfm.loc[dfm['Yfi_val'],['ep_pointing_id','field_id']]

    pointing_group = df_schedule.groupby('ep_pointing_id')

    logging.info(f"EP observations allocated {time_allowed.to(u.hour)}")
    logging.info(f"EP-ZTF schedule covers {len(np.unique(df_schedule['ep_pointing_id']))} EP pointing with {len(df_schedule)} ZTF pointings and {len(np.unique(df_schedule['field_id']))} unique ZTF fields.")



    # construct list queues, pass them back upstream, and have ztf_scheduler run
    # validate_list_queue on them before we rebuild it


    def make_list_queue(fields, time_start_mjd, time_end_mjd, pointing_id):

        # solve TSP
        ordered_fields = check_limits_and_solve_TSP(fields,
            Time(time_start_mjd, format='mjd'), Time(time_end_mjd, format='mjd'))

        targets = []
        for field in ordered_fields:
            targets.append({
                        "program_id": 2,
                        "subprogram_name": "Einstein_Probe",
                        "program_pi": "Kasliwal/Graham",
                        "field_id": field,
                        "max_airmass": 2.5,
                        "filter_id": 2,
                        "n_repeats": 1})
            

        target_df = pd.DataFrame(targets)

        date = Time(time_start_mjd, format='mjd').iso.split()[0]

        data = {}

        data["queue_name"] = f"EP_{date}_{pointing_id}"
        data["queue_description"] = data["queue_name"]
        data["queue_manager"] = "list"
        data["validity_window_mjd"] = [time_start_mjd, time_end_mjd]
        data["targets"] = target_df.to_dict(orient='records')

        # make a fake QueueConfiguration
        queue_config = Configuration(None)
        queue_config.config = data

        logging.info(f"Queue {data['queue_name']} created with {len(ordered_fields)} fields: {ordered_fields}")

        return ListQueueManager(data["queue_name"], queue_config)

    list_queues = []

    for ep_id, grp in pointing_group:

        start_mjd = df_ep_fov.loc[ep_id, 'start_mjd']
        end_mjd = df_ep_fov.loc[ep_id, 'end_mjd']

        unique_fields = np.unique(grp['field_id']).tolist()


        lq = make_list_queue(unique_fields, start_mjd, end_mjd, ep_id) 

        list_queues.append(lq)

    return list_queues
