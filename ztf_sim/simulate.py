from __future__ import print_function
from __future__ import absolute_import
import configparser
import numpy as np
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from .TelescopeStateMachine import TelescopeStateMachine
from .Scheduler import Scheduler
from .QueueManager import GreedyQueueManager, QueueEmptyError, GurobiQueueManager
from .QueueManager import calc_pool_stats, calc_queue_stats
from .ObsLogger import ObsLogger
from .configuration import SchedulerConfiguration
from .constants import BASE_DIR, P48_loc
from .utils import block_index


# check aggressively for setting with copy
import pandas as pd
pd.options.mode.chained_assignment = 'raise'  # default='warn'

# TODO: tag database with commit hash


def simulate(scheduler_config_file, run_config_file = 'default.cfg',
        profile=False, raise_queue_empty=False, fallback=True):

    if profile:
        try:
            from pyinstrument import Profiler
        except ImportError:
            print('Error importing pyinstrument')
            profile = False

    run_config = configparser.ConfigParser()
    run_config_file_fullpath = BASE_DIR+'../config/{}'.format(run_config_file)
    run_config.read(run_config_file_fullpath)

    # load config parameters into local variables
    start_time = run_config['simulation']['start_time']
    try:
        weather_year = run_config['simulation']['weather_year']
    except KeyError:
        weather_year = None
    if (weather_year.lower() == "none"):
        weather_year = None
    else:
        weather_year = int(weather_year)
    survey_duration = \
        run_config['simulation'].getfloat('survey_duration_days') * u.day

    # set up Scheduler
    scheduler_config_file_fullpath = \
            BASE_DIR + '../sims/{}'.format(scheduler_config_file)
    scheduler = Scheduler(scheduler_config_file_fullpath,
            run_config_file_fullpath)
    run_name = scheduler.scheduler_config.config['run_name']

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

    # initialize to a low value so we start by assigning nightly requests
    current_night_mjd = 0

    while tel.current_time < (survey_start_time + survey_duration):

        # check if it is a new night and reload queue with new requests
        if np.floor(tel.current_time.mjd) > current_night_mjd:
            # use the state machine to allow us to skip weathered out nights
            #if tel.check_if_ready():
            scheduler.obs_log.prev_obs = None


            # Look for timed queues that will be valid tonight,
            # to exclude from the nightly solution
            scheduler.timed_queues_tonight = []
            block_start = block_index(tel.current_time)
            block_stop = block_index(Time(tel.current_time.mjd+ 1, format='mjd'))
            exclude_blocks = []
            for qq_name, qq in scheduler.queues.items():
                if qq.queue_name in ['default', 'fallback']:
                    continue
                if qq.validity_window is not None:
                    valid_blocks = qq.valid_blocks(complete_only=True)
                    valid_blocks_tonight = [b for b in valid_blocks if
                            (block_start <= b <= block_stop)]
                    if len(valid_blocks_tonight):
                        scheduler.timed_queues_tonight.append(qq_name)
                    exclude_blocks.extend(valid_blocks_tonight)

            scheduler.Q.assign_nightly_requests(tel.current_state_dict(),
                    scheduler.obs_log, exclude_blocks = exclude_blocks)
            current_night_mjd = np.floor(tel.current_time.mjd)
            # log pool stats
            tel.logger.info(calc_pool_stats(
                scheduler.Q.rp.pool, intro="Nightly requests initialized"))

        if tel.check_if_ready():
            current_state = tel.current_state_dict()

            # drop out of a timed queue if it's no longer valid
            if scheduler.Q.queue_name != 'default':
                if not scheduler.Q.is_valid(current_state['current_time']):
                    scheduler.set_queue('default')


            # check if a timed queue is now valid
            for qq_name in scheduler.timed_queues_tonight:
                qq = scheduler.queues[qq_name]
                if qq.is_valid(current_state['current_time']):
                    # only switch if we are in the default or fallback queue
                    if scheduler.Q.queue_name in ['default', 'fallback']:
                        scheduler.set_queue(qq_name)

            1/0 

            # get coords
            try:
                next_obs = scheduler.Q.next_obs(current_state, 
                        scheduler.obs_log)
                # TODO: debugging check...
                assert(next_obs['request_id'] in scheduler.Q.queue.index)
            except QueueEmptyError:
                if scheduler.Q.queue_name != 'default':
                    scheduler.set_queue('default')
                    try:
                        next_obs = scheduler.Q.next_obs(current_state, 
                                scheduler.obs_log)
                        assert(next_obs['request_id'] in scheduler.Q.queue.index)
                        continue
                    except QueueEmptyError:
                        tel.logger.info("Queue empty!  Trying fallback queue...")
                        if 'fallback' in scheduler.queues:
                            next_obs = scheduler.queues['fallback'].next_obs(
                                    current_state, scheduler.obs_log)
                        else:
                            tel.logger.info("No fallback queue defined!")

                else:
                    if not raise_queue_empty:
                            tel.logger.info("Queue empty!  Waiting...")
                            scheduler.obs_log.prev_obs = None
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
                    scheduler.obs_log.prev_obs = None
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
                scheduler.obs_log.prev_obs = None
                tel.wait()
                continue

            # try to expose
            if not tel.start_exposing():
                tel.set_cant_observe()
                tel.logger.info("Exposure failure!  Waiting...")
                scheduler.obs_log.prev_obs = None
                tel.wait()
                continue
            else:
                # exposure completed successfully.  now
                # a) store exposure information in pointing history sqlite db
                current_state = tel.current_state_dict()
                scheduler.obs_log.log_pointing(current_state, next_obs)
                # b) remove completed request_id from the pool and the queue
                # TODO: debugging check
                assert(next_obs['request_id'] in scheduler.queues[next_obs['queue_name']].queue.index)
                # TODO: check this with request sets...
                scheduler.queues[next_obs['queue_name']].remove_requests(next_obs['request_id']) 
        else:
            scheduler.obs_log.prev_obs = None
            tel.set_cant_observe()
            tel.wait()

    if profile:
        profiler.stop()
        print(profiler.output_text(str=True, color=True))
        with open('../sims/{}_profile.txt'.format(run_name), 'w') as f:
            f.write(profiler.output_text())

    # TODO: gzip logfile
