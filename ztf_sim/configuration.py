"""Classes for parsing scheduler configuration files."""

import pathlib
import json
import numpy as np
import astropy.units as u
from .ObservingProgram import ObservingProgram
from .Fields import Fields
from .constants import PROGRAM_NAMES, PROGRAM_NAME_TO_ID, EXPOSURE_TIME, TIME_BLOCK_SIZE
from .QueueManager import GreedyQueueManager, QueueEmptyError, GurobiQueueManager, ListQueueManager


class Configuration(object):

    def __init__(self, config_file):

        if config_file is not None:
            self.load_configuration(config_file)

    def load_configuration(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
        self.config = config

class SchedulerConfiguration(Configuration):

    def __init__(self, config_file):
        super().__init__(config_file)
        self.scheduler_config_file = pathlib.PurePosixPath(config_file)
        self.check_configuration()

    def check_configuration(self):
        if 'queues' not in self.config:
            raise ValueError("Scheduler configuration must give queues")
        has_default = False
        for queue_pars in self.config['queues']:
            assert "queue_name" in queue_pars
            if queue_pars["queue_name"] == "default":
                has_default = True
            assert "config_file" in queue_pars
        if not has_default:
            raise ValueError("Scheduler configuration must specify a default queue")

    def build_queue_configs(self):

        queue_configs = {}

        for queue_pars in self.config['queues']:
            queue_config = QueueConfiguration(
                    self.scheduler_config_file.parent / queue_pars["config_file"])
            queue_configs[queue_pars["queue_name"]] = queue_config

        return queue_configs

    def build_queues(self, queue_configs):
        
        queues = {}
        for queue_name, queue_config in queue_configs.items():
            
            queue_manager = queue_config.config['queue_manager']
            assert (queue_manager in ('list', 'greedy', 'gurobi'))

            if queue_manager == 'list':
                queues[queue_name] = ListQueueManager(queue_name, queue_config)
            elif queue_manager == 'greedy':
                queues[queue_name] = GreedyQueueManager(queue_name, queue_config)
            elif queue_manager == 'gurobi':
                queues[queue_name] = GurobiQueueManager(queue_name, queue_config)


        return queues





class QueueConfiguration(Configuration):

    def __init__(self, config_file):
        super().__init__(config_file)
        self.check_configuration()

    def check_configuration(self):

        if self.config['queue_manager'] != 'list':
            for month in range(1,13):
                if not np.isclose(np.sum(
                    [prog['program_observing_fraction']*prog['subprogram_fraction'] 
                    for prog in self.config['observing_programs']
                    if ((prog['active_months'] == 'all') or 
                    (month in np.atleast_1d(prog['active_months'])))
                    ]), 1.0):
                    raise ValueError('Observing fractions must sum to 1')

            # could do this via schema validation
            for prog in self.config['observing_programs']:
                if prog['program_name'] not in PROGRAM_NAMES:
                    raise ValueError('{} not in known programs'.format(
                        prog['program_name']))

    def build_observing_programs(self):

        OPs = []
        f = Fields()
        for prog in self.config['observing_programs']:
            assert(('field_ids' in prog) or ('field_selections' in prog))
            assert(('field_ids' in prog) != ('field_selections' in prog))
            if 'field_ids' in prog:
                field_ids = prog['field_ids']
                for field_id in field_ids:
                    if field_id not in f.fields.index:
                        raise ValueError(f'Input field_id {field_id} is not valid')
            else: 
                field_ids = f.select_field_ids(**prog['field_selections'])
            if 'nobs_range' not in prog:
                prog['nobs_range'] = None
            if 'intranight_gap_min' not in prog:
                prog['intranight_gap_min'] = TIME_BLOCK_SIZE
            else:
                # make it a quantity
                prog['intranight_gap_min'] = prog['intranight_gap_min'] * u.minute 
            if 'exposure_time' not in prog:
                prog['exposure_time'] = EXPOSURE_TIME
            else:
                # make it a quantity
                prog['exposure_time'] = prog['exposure_time'] * u.second
            OP = ObservingProgram(PROGRAM_NAME_TO_ID[prog['program_name']],
                                  prog['subprogram_name'], 
                                  prog['program_pi'], 
                                  prog['program_observing_fraction'],
                                  prog['subprogram_fraction'],
                                  field_ids, prog['filter_ids'],
                                  prog['internight_gap_days'] * u.day,
                                  prog['intranight_gap_min'],
                                  prog['n_visits_per_night'],
                                  exposure_time = prog['exposure_time'],
                                  nobs_range = prog['nobs_range'],
                                  filter_choice=prog['filter_choice'],
                                  active_months=prog['active_months'])
            OPs.append(OP)

        return OPs
