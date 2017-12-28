
from builtins import object
import configparser
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .ObsLogger import ObsLogger
from .configuration import ObservingProgramConfiguration
from .constants import BASE_DIR





class Scheduler(object):

    def __init__(self, observing_program_config_file_fullpath, 
            run_config_file_fullpath):

        self.op_config = ObservingProgramConfiguration(
            observing_program_config_file_fullpath)
        self.observing_programs = self.op_config.build_observing_programs()

        self.run_config = configparser.ConfigParser()
        self.run_config.read(run_config_file_fullpath)


        self.Q = None
        # use for swapping queues in the night, e.g. for TOOs
        self.other_queues = {}

        self.set_queue_manager(queue_name = 'default',
            queue_manager = self.run_config['scheduler']['queue_manager'])

        for op in self.observing_programs:
            self.Q.add_observing_program(op)

        

        if 'log_name' in self.run_config['scheduler']:
            log_name = self.run_config['scheduler']['log_name']
        else:
            log_name = self.op_config.config['run_name']

        # initialize sqlite history
        self.obs_log = ObsLogger(log_name,
                clobber=self.run_config['scheduler'].getboolean('clobber_db')) 


    def set_queue_manager(self, queue_name = 'default', 
            queue_manager = 'gurobi', clobber=False):

        # TODO: log the switch

        assert (queue_manager in ('list', 'greedy', 'gurobi'))

        # store previous queue
        if self.Q is not None:
            self.other_queues[self.Q.queue_name] = self.Q

        if clobber or (queue_name not in self.other_queues):
            if queue_manager == 'list':
                self.Q = ListQueueManager(queue_name=queue_name)
            elif queue_manager == 'greedy':
                self.Q = GreedyQueueManager(queue_name=queue_name)
            elif queue_manager == 'gurobi':
                self.Q = GurobiQueueManager(queue_name=queue_name)
        else:
            self.Q = self.other_queues.pop(queue_name)
