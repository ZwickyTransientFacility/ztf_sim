from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from .TelescopeStateMachine import TelescopeStateMachine
from .Scheduler import Scheduler
from .QueueManager import GreedyQueueManager, QueueEmptyError, GurobiQueueManager
from .QueueManager import calc_pool_stats, calc_queue_stats
from .ObsLogger import ObsLogger
from .configuration import ObservingProgramConfiguration
from .constants import BASE_DIR, P48_loc


# check aggressively for setting with copy
import pandas as pd
pd.options.mode.chained_assignment = 'raise'  # default='warn'

# TODO: tag database with commit hash


def simulate(config_file, profile=False, raise_queue_empty=True):

    if profile:
        try:
            from pyinstrument import Profiler
        except ImportError:
            print('Error importing pyinstrument')
            profile = False

    ztf_config = ObservingProgramConfiguration(
            BASE_DIR + '../sims/{}'.format(config_file))

    # load config parameters into local variables
    run_name = ztf_config.config['run_name']
    start_time = ztf_config.config['start_time']
    weather_year = ztf_config.config['weather_year']
    if weather_year == "None":
        weather_year = None
    survey_duration = ztf_config.config['survey_duration_days'] * u.day
    block_programs = ztf_config.config['block_programs']
    observing_programs = ztf_config.build_observing_programs()

    if profile:
        if survey_duration > 1. * u.day:
            print("Don't profile long runs: 25% overhead")
            profile = False
        else:
            profiler = Profiler()
            profiler.start()

    survey_start_time = Time(start_time, scale='utc', location=P48_loc)

    tel = TelescopeStateMachine(
        current_time=survey_start_time,
        historical_observability_year=weather_year,
        logfile=BASE_DIR + '../sims/{}_log.txt'.format(run_name))

    # set up Scheduler
    # TODO: make queue manager a config parameter
    scheduler = Scheduler(config_file, queue_manager='gurobi')

    # initialize nightly field requests (Tom Barlow function)
    scheduler.Q.assign_nightly_requests(tel.current_state_dict())
    # log pool stats
    tel.logger.info(calc_pool_stats(
        scheduler.Q.rp.pool, intro="Nightly requests initialized"))

    # initialize sqlite history
    log = ObsLogger(run_name, tel.current_time)

    current_night_mjd = np.floor(tel.current_time.mjd)

    while tel.current_time < (survey_start_time + survey_duration):

        # check if it is a new night and reload queue with new requests
        if np.floor(tel.current_time.mjd) > current_night_mjd:
            # use the state machine to allow us to skip weathered out nights
            log.prev_obs = None
            scheduler.Q.assign_nightly_requests(tel.current_state_dict())
            current_night_mjd = np.floor(tel.current_time.mjd)
            # log pool stats
            tel.logger.info(calc_pool_stats(
                scheduler.Q.rp.pool, intro="Nightly requests initialized"))

        if tel.check_if_ready():
            current_state = tel.current_state_dict()
            # get coords
            try:
                next_obs = scheduler.Q.next_obs(current_state)
                # TODO: debugging check...
                assert(next_obs['request_id'] in scheduler.Q.queue.index)
            except QueueEmptyError:
                if not raise_queue_empty:
                    tel.logger.info("Queue empty!  Waiting...")
                    log.prev_obs = None
                    tel.wait()
                    continue
                else:
                    tel.logger.info(calc_queue_stats(
                        scheduler.Q.queue, current_state,
                        intro="Queue returned no next_obs. Current queue status:"))
                    tel.logger.info(calc_pool_stats(
                        scheduler.Q.rp.pool, intro="Current pool status:"))

                    # TODO: in py3, chained exceptions come for free
                    raise QueueEmptyError

            # try to change filters, if needed
            if next_obs['target_filter_id'] != current_state['current_filter_id']:
                if not tel.start_filter_change(next_obs['target_filter_id']):
                    tel.logger.info("Filter change failure!  Waiting...")
                    log.prev_obs = None
                    tel.wait()
                    continue

            # try to slew to the next target
            if not tel.start_slew(coord.SkyCoord(next_obs['target_ra'] * u.deg,
                                                 next_obs['target_dec'] * u.deg)):
                tel.set_cant_observe()
                # TODO: log the failure
                # "missed history": http://ops2.lsst.org/docs/current/architecture.html#output-tables
                tel.logger.info("Failure slewing to {}, {}!  Waiting...".format
                                (next_obs['target_ra'] * u.deg, next_obs['target_dec'] * u.deg))
                log.prev_obs = None
                tel.wait()
                continue

            # try to expose
            if not tel.start_exposing():
                tel.set_cant_observe()
                tel.logger.info("Exposure failure!  Waiting...")
                log.prev_obs = None
                tel.wait()
                continue
            else:
                # exposure completed successfully.  now
                # a) store exposure information in pointing history sqlite db
                current_state = tel.current_state_dict()
                log.log_pointing(current_state, next_obs)
                # b) update Fields
                scheduler.Q.fields.mark_field_observed(next_obs, 
                        current_state['current_time'])
                # c) remove completed request_id from the pool and the queue
                # TODO: debugging check
                assert(next_obs['request_id'] in scheduler.Q.queue.index)
                # TODO: check this with request sets...
                scheduler.Q.remove_requests(next_obs['request_id'])
        else:
            log.prev_obs = None
            tel.set_cant_observe()
            tel.wait()

    if profile:
        profiler.stop()
        print(profiler.output_text(str=True, color=True))
        with open('../sims/{}_profile.txt'.format(run_name), 'w') as f:
            f.write(profiler.output_text())

    # TODO: gzip logfile