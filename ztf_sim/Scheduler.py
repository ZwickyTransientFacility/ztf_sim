
from builtins import object
import configparser
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .ObsLogger import ObsLogger
from .configuration import SchedulerConfiguration
from .constants import BASE_DIR





class Scheduler(object):

    def __init__(self, scheduler_config_file_fullpath, 
            run_config_file_fullpath, other_queue_configs = None):

        self.scheduler_config = SchedulerConfiguration(
            scheduler_config_file_fullpath)
        self.queue_configs = self.scheduler_config.build_queue_configs()
        self.queues = self.scheduler_config.build_queues(self.queue_configs)

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
            del self.queues[queue_name] 
        else:
            raise ValueError(f"Queue {queue_name} does not exist!")
