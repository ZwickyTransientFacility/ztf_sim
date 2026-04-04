"""Classes for parsing scheduler configuration files."""

import pathlib
import json
import numpy as np
import astropy.units as u
from .ObservingProgram import ObservingProgram
from .Fields import Fields
from .constants import PROGRAM_NAMES, PROGRAM_NAME_TO_ID, EXPOSURE_TIME, TIME_BLOCK_SIZE
from .QueueManager import GreedyQueueManager, QueueEmptyError, GurobiQueueManager, ListQueueManager
from .field_selection_functions import *


class Configuration(object):
    """Base class for loading JSON scheduler configuration files.

    Attributes
    ----------
    config : dict
        Parsed JSON configuration.
    """

    def __init__(self, config_file):
        """Load a JSON configuration file.

        Parameters
        ----------
        config_file : str or pathlib.Path or None
            Path to the JSON configuration file. If ``None``, no file is
            loaded (used by `MMASkymap` to build synthetic configurations).
        """

        if config_file is not None:
            self.load_configuration(config_file)

    def load_configuration(self, config_file):
        """Parse a JSON file and store the result in ``self.config``.

        Parameters
        ----------
        config_file : str or pathlib.Path
            Path to the JSON configuration file.
        """
        with open(config_file, 'r') as f:
            config = json.load(f)
        self.config = config

class SchedulerConfiguration(Configuration):
    """Top-level scheduler configuration specifying which queues to run.

    Parses a JSON file that lists one or more queue configurations. Exactly
    one queue must be named ``'default'``.
    """

    def __init__(self, config_file):
        """Load and validate the scheduler configuration.

        Parameters
        ----------
        config_file : str or pathlib.Path
            Path to the JSON scheduler configuration file.

        Raises
        ------
        ValueError
            If ``'queues'`` is absent or no queue is named ``'default'``.
        """
        super().__init__(config_file)
        self.scheduler_config_file = pathlib.PurePosixPath(config_file)
        self.check_configuration()

    def check_configuration(self):
        """Validate the scheduler configuration.

        Raises
        ------
        ValueError
            If the ``'queues'`` key is missing or no queue entry has
            ``queue_name == 'default'``.
        AssertionError
            If any queue entry is missing ``'queue_name'`` or
            ``'config_file'``.
        """
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
        """Load each queue's configuration file into a `QueueConfiguration`.

        Returns
        -------
        dict
            Mapping ``queue_name (str) -> QueueConfiguration``.

        Raises
        ------
        Exception
            Re-raised from `QueueConfiguration` if a config file cannot be
            read or parsed.
        """

        queue_configs = {}

        for queue_pars in self.config['queues']:
            try:
                queue_config = QueueConfiguration(
                    self.scheduler_config_file.parent / queue_pars["config_file"])
            except Exception as e:
                print(f'Error reading config file {queue_pars["config_file"]}')
                raise(e)

            queue_configs[queue_pars["queue_name"]] = queue_config

        return queue_configs

    def build_queues(self, queue_configs):
        """Instantiate queue managers from their configurations.

        Parameters
        ----------
        queue_configs : dict
            Mapping ``queue_name -> QueueConfiguration`` as returned by
            `build_queue_configs`.

        Returns
        -------
        dict
            Mapping ``queue_name (str) -> QueueManager`` subclass instance
            (``ListQueueManager``, ``GreedyQueueManager``, or
            ``GurobiQueueManager``).

        Raises
        ------
        AssertionError
            If ``queue_manager`` value is not one of ``'list'``,
            ``'greedy'``, ``'gurobi'``.
        """

        queues = {}
        for queue_name, queue_config in queue_configs.items():
            
            queue_manager = queue_config.config['queue_manager']
            assert (queue_manager in ('list', 'greedy', 'gurobi'))

            try:
                if queue_manager == 'list':
                    queues[queue_name] = ListQueueManager(queue_name, queue_config)
                elif queue_manager == 'greedy':
                    queues[queue_name] = GreedyQueueManager(queue_name, queue_config)
                elif queue_manager == 'gurobi':
                    queues[queue_name] = GurobiQueueManager(queue_name, queue_config)
            except Exception as e:
                print(f'Error building queue {queue_name}')
                raise(e)


        return queues





class QueueConfiguration(Configuration):
    """Per-queue configuration specifying observing programs and scheduling mode.

    Parses the JSON file for a single queue. Validates that observing fractions
    sum to 1 for every active month and that all program names are recognised.
    """

    def __init__(self, config_file):
        """Load and validate a queue configuration file.

        Parameters
        ----------
        config_file : str or pathlib.Path
            Path to the JSON queue configuration file.

        Raises
        ------
        ValueError
            If observing fractions do not sum to 1 for any active month, or
            if an unknown program name is encountered.
        """
        super().__init__(config_file)
        self.check_configuration()

    def check_configuration(self):
        """Validate observing fractions and program names.

        Raises
        ------
        ValueError
            If ``program_observing_fraction × subprogram_fraction`` values do
            not sum to 1.0 for any calendar month in which at least one
            program is active, or if a program name is not in
            ``PROGRAM_NAMES``.
        """

        if self.config['queue_manager'] != 'list' and len(self.config['observing_programs']):
            for month in range(1,13):
                op_sum = np.sum(
                    [prog['program_observing_fraction']*prog['subprogram_fraction'] 
                    for prog in self.config['observing_programs']
                    if ((prog['active_months'] == 'all') or 
                    (month in np.atleast_1d(prog['active_months'])))
                    ])
                if not np.isclose(op_sum, 1.0):
                    raise ValueError(f"Observing fractions must sum to 1 ({op_sum}): {[(prog['subprogram_name'], prog['program_observing_fraction']*prog['subprogram_fraction']) for prog in self.config['observing_programs']]}")

            # could do this via schema validation
            for prog in self.config['observing_programs']:
                if prog['program_name'] not in PROGRAM_NAMES:
                    raise ValueError('{} not in known programs'.format(
                        prog['program_name']))

    def build_observing_programs(self):
        """Instantiate `ObservingProgram` objects from the queue configuration.

        Resolves field sources: converts ``field_selections`` dicts to
        explicit field ID lists via `Fields.select_field_ids`, validates
        ``field_ids`` against the field grid, and checks that
        ``field_selection_function`` names exist in
        ``field_selection_functions``.

        Returns
        -------
        list of ObservingProgram
            One entry per program defined in ``self.config['observing_programs']``.

        Raises
        ------
        ValueError
            If a ``field_ids`` entry is not a valid ZTF field ID.
        AssertionError
            If a program provides more than one field source, or if a
            ``field_selection_function`` name is not defined.
        """

        OPs = []
        f = Fields()
        for prog in self.config['observing_programs']:
            assert(('field_ids' in prog) or ('field_selections' in prog)
                    or ('field_selection_function' in prog))
            assert(('field_ids' in prog) + ('field_selections' in prog) + 
                    ('field_selection_function' in prog) == 1)
            if 'field_ids' in prog:
                field_ids = prog['field_ids']
                for field_id in field_ids:
                    if field_id not in f.fields.index:
                        raise ValueError(f'Input field_id {field_id} is not valid')
                field_selection_function = None
            elif 'field_selections' in prog: 
                field_ids = f.select_field_ids(**prog['field_selections'])
                field_selection_function = None
            else:
                field_ids = None
                field_selection_function = prog['field_selection_function']
                # check if it exists
                # special case the EP selection, which is done by 
                # make_nightly_timed_blocks
                if field_selection_function != 'EP-bypass':
                    assert(field_selection_function in globals())
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
                                  active_months=prog['active_months'],
                                  field_selection_function = field_selection_function)
            OPs.append(OP)

        return OPs
