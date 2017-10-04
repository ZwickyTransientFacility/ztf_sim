
from builtins import object
import configparser
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .ObsLogger import ObsLogger
from .configuration import ObservingProgramConfiguration
from .constants import BASE_DIR





class Scheduler(object):

    def __init__(self, observing_program_config_file, run_config_file, 
            queue_manager = 'gurobi'):

        self.op_config = ObservingProgramConfiguration(
            BASE_DIR + '../sims/{}'.format(observing_program_config_file))
        self.observing_programs = self.op_config.build_observing_programs()

        self.run_config = configparser.ConfigParser()
        self.run_config.read(BASE_DIR + '../config/{}'.format(run_config_file))


        self.Q = None
        # use for swapping queues in the night, e.g. for TOOs
        self.prev_Q = None

        self.set_queue_manager(
            queue_manager = self.run_config['scheduler']['queue_manager'])

        for op in self.observing_programs:
            self.Q.add_observing_program(op)

        

        if 'log_name' in self.run_config['scheduler']:
            log_name = self.run_config['scheduler']['log_name']
        else:
            log_name = self.op_config['run_name']

        # initialize sqlite history
        self.log = ObsLogger(log_name,
                clobber=self.run_config['scheduler'].getboolean('clobber_db')) 


    def set_queue_manager(self, queue_manager = 'gurobi'):

        assert (queue_manager in ('list', 'greedy', 'gurobi'))

        self.prev_Q = self.Q

        if queue_manager == 'list':
            self.Q = ListQueueManager()
        elif queue_manager == 'greedy':
            self.Q = GreedyQueueManager()
        elif queue_manager == 'gurobi':
            self.Q = GurobiQueueManager()

    def pop_prev_queue(self):

        self.Q = self.prev_Q
