
from builtins import object
import configparser
import numpy as np
from astropy.time import Time
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .ObsLogger import ObsLogger
from .configuration import SchedulerConfiguration
from .constants import BASE_DIR
from .utils import block_index





class Scheduler(object):

    def __init__(self, scheduler_config_file_fullpath, 
            run_config_file_fullpath, other_queue_configs = None):

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
                clobber=self.run_config['scheduler'].getboolean('clobber_db')) 


    def set_queue(self, queue_name): 

        # TODO: log the switch
        
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

    def find_excluded_blocks_tonight(self, time_now):
        # also sets up timed_queues_tonight

        # start of the night
        mjd_today = np.floor(time_now.mjd).astype(int)

        # Look for timed queues that will be valid tonight,
        # to exclude from the nightly solution
        self.timed_queues_tonight = []
        block_start = block_index(Time(mjd_today, format='mjd'))
        block_stop = block_index(Time(mjd_today + 1, format='mjd'))
        exclude_blocks = []
        for qq_name, qq in self.queues.items():
            if qq.queue_name in ['default', 'fallback']:
                continue
            if qq.validity_window is not None:
                valid_blocks = qq.valid_blocks(complete_only=True)
                all_valid_blocks = qq.valid_blocks(complete_only=False)
                valid_blocks_tonight = [b for b in valid_blocks if
                        (block_start <= b <= block_stop)]
                all_valid_blocks_tonight = [b for b in all_valid_blocks if
                        (block_start <= b <= block_stop)]
                if len(all_valid_blocks_tonight):
                    self.timed_queues_tonight.append(qq_name)
                exclude_blocks.extend(valid_blocks_tonight)
        return exclude_blocks

    def check_for_TOO_queue_and_switch(self, time_now):
        # check if a TOO queue is now valid
        for qq_name, qq in self.queues.items():
            if qq.is_TOO:
                if qq.is_valid(time_now):
                    # only switch if we don't have an active TOO queue
                    if not self.Q.is_TOO and len(qq.queue):
                        self.set_queue(qq_name)

    def check_for_timed_queue_and_switch(self, time_now):
            # drop out of a timed queue if it's no longer valid
            if self.Q.queue_name != 'default':
                if not self.Q.is_valid(time_now):
                    self.set_queue('default')

            # check if a timed queue is now valid
            for qq_name in self.timed_queues_tonight:
                qq = self.queues[qq_name]
                if qq.is_valid(time_now): 
                    if (qq.queue_type == 'list'): 
                        # list queues should have items in them
                        if len(qq.queue):
                            # only switch if we are in the default or fallback queue
                            if self.Q.queue_name in ['default', 'fallback']:
                                self.set_queue(qq_name)
                    else:
                        # don't have a good way to check length of non-list
                        # queues before nightly assignments
                        if self.Q.queue_name in ['default', 'fallback']:
                            self.set_queue(qq_name)

    def remove_empty_and_expired_queues(self, time_now):
        queues_for_deletion = []
        for qq_name, qq in self.queues.items():
            if qq.queue_name in ['default', 'fallback']:
                continue
            if qq.validity_window is not None:
                if qq.validity_window[1] < time_now:
                    queues_for_deletion.append(qq_name)
                    continue
            if len(qq.queue) == 0:
                    queues_for_deletion.append(qq_name)

        # ensure we don't have duplicate values
        queues_for_deletion = set(queues_for_deletion)

        for qq_name in queues_for_deletion:
            self.delete_queue(qq_name)
