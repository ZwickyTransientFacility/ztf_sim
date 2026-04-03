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
from .field_selection.ep import make_ep_blocks





class Scheduler(object):
    """Top-level scheduler that owns all queues, the observation log, and skymaps.

    Orchestrates nightly queue assignment, Einstein Probe timed-block setup,
    queue switching (TOO and timed), and accounting for timed program time
    commitments.

    Attributes
    ----------
    queues : dict
        Mapping ``queue_name -> QueueManager`` for all active queues.
    Q : QueueManager
        Currently active queue.
    obs_log : ObsLogger
        Observation history and output database.
    skymaps : dict
        Mapping ``trigger_name -> MMASkymap`` for registered skymaps.
    timed_queues_tonight : list of str
        Queue names that have valid timed windows tonight.
    """

    def __init__(self, scheduler_config_file_fullpath,
            run_config_file_fullpath, other_queue_configs=None,
            output_path=BASE_DIR+'../sims/'):
        """Initialise the scheduler from configuration files.

        Parameters
        ----------
        scheduler_config_file_fullpath : str
            Absolute path to the JSON scheduler configuration file.
        run_config_file_fullpath : str
            Absolute path to the INI simulation/run configuration file.
            Must contain a ``[scheduler]`` section with at least
            ``clobber_db``; optionally ``log_name``.
        other_queue_configs : dict or None, optional
            Additional ``{queue_name: QueueConfiguration}`` entries to merge
            into the built queues. Reserved for future use. Default is
            ``None``.
        output_path : str, optional
            Directory where the output SQLite database and log file are
            written. Default is ``../sims/`` relative to the package root.
        """

        self.logger = logging.getLogger(__name__)

        self.scheduler_config = SchedulerConfiguration(
            scheduler_config_file_fullpath)
        self.queue_configs = self.scheduler_config.build_queue_configs()
        self.queues = self.scheduler_config.build_queues(self.queue_configs)
        self.timed_queues_tonight = []

        # used to trigger nightly recomputes
        self.mjd_today = 0

        self.skymaps = {}

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

    def assign_nightly_requests(self, current_state_dict,
                                time_limit=15.*u.minute):
        """Trigger the nightly ILP scheduling for the default queue.

        Computes block-level time commitments from timed queues and timed
        observation counts, then delegates to the default queue's
        ``assign_nightly_requests``.

        Parameters
        ----------
        current_state_dict : dict
            Telescope state dict as returned by
            ``TelescopeStateMachine.current_state_dict()``.
        time_limit : astropy.units.Quantity, optional
            Gurobi wall-clock time limit for the ILP solver. Default is
            15 minutes.
        """
        # Look for timed queues that will be valid tonight,
        # to exclude from the nightly solution
        block_use = self.find_block_use_tonight(current_state_dict['current_time'])
        timed_obs_count = self.count_timed_observations_tonight()

        self.logger.info(f'Block use by timed queues: {block_use}')

        self.queues['default'].assign_nightly_requests(
                        current_state_dict,
                        self.obs_log, 
                        time_limit = time_limit,
                        block_use = block_use,
                        timed_obs_count = timed_obs_count,
                        skymaps = self.skymaps)

    def make_nightly_timed_blocks(self, current_state_dict,
                                time_limit=15.*u.minute):
        """Set up time-windowed queues for tonight, currently Einstein Probe only.

        Clears any EP queues left over from the previous night, checks whether
        the Einstein Probe programme has a non-zero allocation tonight, and if
        so calls `field_selection.ep.make_ep_blocks` to create
        ``ListQueueManager`` objects for each EP observation window. The
        resulting queues are added to ``self.queues`` and tracked in
        ``self.timed_queues_tonight``.

        Parameters
        ----------
        current_state_dict : dict
            Telescope state dict as returned by
            ``TelescopeStateMachine.current_state_dict()``.
        time_limit : astropy.units.Quantity, optional
            Gurobi time limit passed to ``make_ep_blocks``. Default is
            15 minutes.
        """

        to_delete = []
        # clean out any lingering EP queues.
        # don't delete in the loop to avoid "dictionary changed size during iteration" error
        for queue_name in self.queues.keys():
            if queue_name.startswith('EP_20'):
                to_delete.append(queue_name)

        for queue_name in to_delete:
            self.delete_queue(queue_name)

        # this won't (yet) include the EP observations
        block_use = self.find_block_use_tonight(current_state_dict['current_time'])
        timed_obs_count = self.count_timed_observations_tonight()

        # this is duplicative but we need it
        # don't worry about counting timed obs since it will be redone in 
        # assign_nightly_requests
        self.queues['default'].determine_allowed_requests(
                    current_state_dict['current_time'],
                    self.obs_log, timed_obs_count = timed_obs_count)

        programs = self.queues['default'].requests_allowed.keys()

        if (3, 'Einstein_Probe') not in programs:
            logging.info(f'No additional nightly timed queues to add.')
            return
        
        # TODO: make this more general if desired
        time_allowed = self.queues['default'].requests_allowed[
                (3, 'Einstein_Probe')] * (30 * u.second + READOUT_TIME)

        if time_allowed < 0*u.second:
            logging.info(f'No time available for new timed queues.')

        other_timed_queues_tonight = [self.queues[qq] for qq in self.timed_queues_tonight]

        ep_queues = make_ep_blocks(current_state_dict['current_time'],
                                   time_allowed,
                                   time_limit=time_limit,
                                   other_timed_queues_tonight=other_timed_queues_tonight)

        for eq in ep_queues:
            self.add_queue(eq.queue_name, eq)
            logging.info(f'Added queue {eq.queue_name}')

    def set_queue(self, queue_name):
        """Switch the active queue to *queue_name*.

        Parameters
        ----------
        queue_name : str
            Name of the queue to activate. Must exist in ``self.queues``.

        Raises
        ------
        ValueError
            If *queue_name* is not in ``self.queues``.
        """

        if queue_name not in self.queues:
            raise ValueError(f'Requested queue {queue_name} not available!')

        self.Q = self.queues[queue_name]
        
    def add_queue(self, queue_name, queue, clobber=True):
        """Register a new queue with the scheduler.

        Parameters
        ----------
        queue_name : str
            Key under which to store the queue.
        queue : QueueManager
            Queue instance to register.
        clobber : bool, optional
            If ``True`` (default), silently replace any existing queue with
            the same name. If ``False``, raise ``ValueError`` on conflict.

        Raises
        ------
        ValueError
            If *clobber* is ``False`` and *queue_name* already exists.
        """

        if clobber or (queue_name not in self.queues):
            self.queues[queue_name] = queue 
        else:
            raise ValueError(f"Queue {queue_name} already exists!")

    def delete_queue(self, queue_name):
        """Remove a queue from the scheduler.

        If the deleted queue is the currently active queue, the scheduler
        automatically switches to the ``'default'`` queue.

        Parameters
        ----------
        queue_name : str
            Name of the queue to delete.

        Raises
        ------
        ValueError
            If *queue_name* does not exist in ``self.queues``.
        """

        if (queue_name in self.queues):
            if self.Q.queue_name == queue_name:
                self.set_queue('default')
            del self.queues[queue_name] 
        else:
            raise ValueError(f"Queue {queue_name} does not exist!")

    def add_skymap(self, trigger_name, skymap, clobber=True):
        """Register a multi-messenger skymap with the scheduler.

        Parameters
        ----------
        trigger_name : str
            Unique identifier for the event (e.g. ``'S190814bv'``).
        skymap : MMASkymap
            Skymap object to register.
        clobber : bool, optional
            If ``True`` (default), replace any existing skymap with the same
            name. If ``False``, raise ``ValueError`` on conflict.

        Raises
        ------
        ValueError
            If *clobber* is ``False`` and *trigger_name* already exists.
        """

        if clobber or (trigger_name not in self.skymaps):
            self.skymaps[trigger_name] = skymap
        else:
            raise ValueError(f"Skymap {trigger_name} already exists!")

    def delete_skymap(self, trigger_name):
        """Remove a skymap from the scheduler.

        Parameters
        ----------
        trigger_name : str
            Name of the skymap to remove.

        Raises
        ------
        ValueError
            If *trigger_name* does not exist in ``self.skymaps``.
        """

        if (trigger_name in self.skymaps):
            del self.skymaps[trigger_name] 
        else:
            raise ValueError(f"Skymap {trigger_name} does not exist!")

    def find_block_use_tonight(self, time_now):
        """Compute the fraction of each block already committed to timed queues.

        Also populates ``self.timed_queues_tonight`` with the names of queues
        that have valid windows tonight.

        Parameters
        ----------
        time_now : astropy.time.Time
            Current simulation time (used to identify tonight).

        Returns
        -------
        block_use : collections.defaultdict
            Mapping ``block_index -> fraction_used`` (0–1). Includes the
            fractions of the twilight blocks consumed by daytime, plus any
            fractions consumed by timed queues.
        """
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
        """Count equivalent standard exposures in timed queues for each programme.

        Returns
        -------
        dict
            Mapping ``program_id -> int`` equivalent-observation count summed
            across all timed queues valid tonight. Programmes with no timed
            observations return 0.
        """
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
        """Switch to an active Target-of-Opportunity queue if one is available.

        Activates the first valid TOO queue with a non-empty observation list.
        If the current queue is already a TOO queue and it has become empty,
        switches to the next available TOO queue.

        Parameters
        ----------
        time_now : astropy.time.Time
            Current simulation time.
        """
        # check if a TOO queue is now valid
        for qq_name, qq in self.queues.items():
            if qq.is_TOO:
                if qq.is_valid(time_now):
                    # switch if the current queue is not a TOO
                    if (not self.Q.is_TOO) and len(qq.queue):
                        self.set_queue(qq_name)
                    # or if the current TOO queue is empty
                    if ((self.Q.is_TOO) and (len(self.Q.queue) == 0) 
                            and len(qq.queue)):
                        self.set_queue(qq_name)

    def check_for_timed_queue_and_switch(self, time_now):
        """Switch into or out of a timed queue based on validity windows.

        Drops back to the ``'default'`` queue if the current non-default queue
        is no longer valid. From the default or fallback queue, switches to
        the first valid timed queue that has observations available.

        Parameters
        ----------
        time_now : astropy.time.Time
            Current simulation time.
        """
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
        """Delete expired or empty non-default queues.

        A queue is removed if its validity window has passed or if it is a
        ``ListQueueManager`` with no remaining observations.

        Parameters
        ----------
        time_now : astropy.time.Time
            Current simulation time used to check validity windows.
        """
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
