
from builtins import object
from .QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from .configuration import ObservingProgramConfiguration
from .constants import BASE_DIR





class Scheduler(object):

    def __init__(self, config_file, queue_manager = 'gurobi'):

        self.config_file = config_file
        self.config = ObservingProgramConfiguration(
            BASE_DIR + '../sims/{}'.format(config_file))
        self.observing_programs = self.config.build_observing_programs()

        self.Q = None
        # use for swapping queues in the night, e.g. for TOOs
        self.prev_Q = None

        self.set_queue_manager(queue_manager = queue_manager)

        for op in self.observing_programs:
            self.Q.add_observing_program(op)


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
