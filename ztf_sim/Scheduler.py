"""Core scheduler classes."""

import configparser
from collections import defaultdict
import logging
import numpy as np
from astropy.time import Time
import astropy.units as u
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .ObsLogger import ObsLogger
from .configuration import SchedulerConfiguration
from .constants import BASE_DIR, PROGRAM_IDS, EXPOSURE_TIME, READOUT_TIME
from .utils import block_index, block_use_fraction
from .utils import next_12deg_evening_twilight, next_12deg_morning_twilight





class Scheduler(object):

    def __init__(self, scheduler_config_file_fullpath, 
            run_config_file_fullpath, other_queue_configs = None,
            output_path = BASE_DIR+'../sims/'):

        self.logger = logging.getLogger(__name__)

        self.scheduler_config = SchedulerConfiguration(
            scheduler_config_file_fullpath)
        self.queue_configs = self.scheduler_config.build_queue_configs()
        self.queues = self.scheduler_config.build_queues(self.queue_configs)
        self.timed_queues_tonight = []

        self.set_queue('default')
        
        self.run_config = configparser.ConfigParser()
        self.run_config.read(run_config_file_fullpath)

        if 'log_name' in self.run_config['scheduler']:
            log_name = self.run_config['scheduler']['log_name']
        else:
            log_name = self.scheduler_config.config['run_name']

        # initialize sqlite history
        self.obs_log = ObsLogger(log_name,
                output_path = output_path,
                clobber=self.run_config['scheduler'].getboolean('clobber_db'),) 


    def set_queue(self, queue_name): 

        if queue_name not in self.queues:
            raise ValueError(f'Requested queue {queue_name} not available!')

        self.Q = self.queues[queue_name]
        


    def add_queue(self,  queue_name, queue, clobber=True):

        if clobber or (queue_name not in self.queues):
            self.queues[queue_name] = queue 
        else:
            raise ValueError(f"Queue {queue_name} already exists!")

    def delete_queue(self, queue_name):

        if (queue_name in self.queues):
            if self.Q.queue_name == queue_name:
                self.set_queue('default')
            del self.queues[queue_name] 
        else:
            raise ValueError(f"Queue {queue_name} does not exist!")

    def find_block_use_tonight(self, time_now):
        # also sets up timed_queues_tonight

        # start of the night
        mjd_today = np.floor(time_now.mjd).astype(int)

        # Look for timed queues that will be valid tonight,
        # to exclude from the nightly solution
        self.timed_queues_tonight = []
        today = Time(mjd_today, format='mjd')
        tomorrow = Time(mjd_today + 1, format='mjd')
        block_start = block_index(today)[0]
        block_stop = block_index(tomorrow)[0]

        block_use = defaultdict(float)

        # compute fraction of twilight blocks not available
        evening_twilight = next_12deg_evening_twilight(today)
        morning_twilight = next_12deg_morning_twilight(today)

        evening_twilight_block = block_index(evening_twilight)[0]
        frac_evening_twilight = block_use_fraction(
                evening_twilight_block, today, evening_twilight)
        block_use[evening_twilight_block] = frac_evening_twilight
        self.logger.debug(f'{frac_evening_twilight} of block {evening_twilight_block} is before 12 degree twilight')

        morning_twilight_block = block_index(morning_twilight)[0]
        frac_morning_twilight = block_use_fraction(
                morning_twilight_block, morning_twilight, tomorrow)
        block_use[morning_twilight_block] = frac_morning_twilight
        self.logger.debug(f'{frac_morning_twilight} of block {morning_twilight_block} is before 12 degree twilight')

        for qq_name, qq in self.queues.items():
            if qq.queue_name in ['default', 'fallback']:
                continue
            if qq.validity_window is not None:
                qq_block_use = qq.compute_block_use()

                is_tonight = False

                # sum block use
                for block, frac in qq_block_use.items():
                    if (block_start <= block <= block_stop):
                        if frac > 0:
                            is_tonight = True
                        self.logger.debug(f'{frac} of block {block} used by queue {qq.queue_name}')
                        block_use[block] += frac
                        if block_use[block] > 1:
                            self.logger.warn(f'Too many observations for block {block}: {block_use[block]}')
                            block_use[block] = 1.

                if is_tonight:    
                    self.timed_queues_tonight.append(qq_name)

        return block_use

    def count_timed_observations_tonight(self):
        # determine how many equivalent obs are in timed queues
        
        timed_obs = {prog:0 for prog in PROGRAM_IDS} 
        if len(self.timed_queues_tonight) == 0:
            return timed_obs

        for qq in self.timed_queues_tonight:
            queue = self.queues[qq].queue.copy()
            if 'n_repeats' not in queue.columns:
                queue['n_repeats'] = 1.
            queue['total_time'] = (queue['exposure_time'] + 
                READOUT_TIME.to(u.second).value)*queue['n_repeats']
            net = queue[['program_id','total_time']].groupby('program_id').agg(np.sum)
            count_equivalent = np.round(net['total_time']/(EXPOSURE_TIME + READOUT_TIME).to(u.second).value).astype(int).to_dict()
            for k, v in count_equivalent.items():
                timed_obs[k] += v

        return timed_obs

    def check_for_TOO_queue_and_switch(self, time_now):
        # check if a TOO queue is now valid
        for qq_name, qq in self.queues.items():
            if qq.is_TOO:
                if qq.is_valid(time_now):
                    # don't switch if there's already a TOO queue active
                    if (not self.Q.is_TOO) and len(qq.queue):
                        self.set_queue(qq_name)

    def check_for_timed_queue_and_switch(self, time_now):
            # drop out of a timed queue if it's no longer valid
            if self.Q.queue_name != 'default':
                if not self.Q.is_valid(time_now):
                    self.set_queue('default')

            # only switch from default or fallback queues
            if self.Q.queue_name in ['default', 'fallback']:
                # check if a timed queue is now valid
                for qq_name, qq in self.queues.items():
                    if (qq.validity_window is not None) and (qq.is_valid(time_now)): 
                        if (qq.queue_type == 'list'): 
                            # list queues should have items in them
                            if len(qq.queue):
                                self.set_queue(qq_name)
                        else:
                            # don't have a good way to check length of non-list
                            # queues before nightly assignments
                            if qq.requests_in_window:
                                self.set_queue(qq_name)

    def remove_empty_and_expired_queues(self, time_now):
        queues_for_deletion = []
        for qq_name, qq in self.queues.items():
            if qq.queue_name in ['default', 'fallback']:
                continue
            if qq.validity_window is not None:
                if qq.validity_window[1] < time_now:
                    self.logger.info(f'Deleting expired queue {qq_name}')
                    queues_for_deletion.append(qq_name)
                    continue
            if (qq.queue_type == 'list') and (len(qq.queue) == 0):
                    self.logger.info(f'Deleting empty queue {qq_name}')
                    queues_for_deletion.append(qq_name)

        # ensure we don't have duplicate values
        queues_for_deletion = set(queues_for_deletion)

        for qq_name in queues_for_deletion:
            self.delete_queue(qq_name)
