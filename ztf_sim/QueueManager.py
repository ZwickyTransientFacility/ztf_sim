"""Queue classes."""

import os
from collections import defaultdict
from datetime import datetime
import logging
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time, TimeDelta
import astroplan
from .Fields import Fields
from .optimize import tsp_optimize, night_optimize
from .cadence import enough_gap_since_last_obs
from .constants import P48_loc, PROGRAM_IDS, FILTER_IDS, TIME_BLOCK_SIZE
from .constants import EXPOSURE_TIME, READOUT_TIME, FILTER_CHANGE_TIME, slew_time
from .constants import PROGRAM_BLOCK_SEQUENCE, LEN_BLOCK_SEQUENCE, MAX_AIRMASS
from .constants import BASE_DIR
from .utils import approx_hours_of_darkness
from .utils import skycoord_to_altaz, seeing_at_pointing
from .utils import altitude_to_airmass, airmass_to_altitude, RA_to_HA, HA_to_RA
from .utils import scalar_len, nightly_blocks, block_index, block_index_to_time
from .utils import block_use_fraction, maximum_altitude, compute_limiting_mag

class QueueEmptyError(Exception):
    """Error class for when the nightly queue has no more fields"""
    pass


class QueueManager(object):

    def __init__(self, queue_name, queue_configuration, rp=None, fields=None):

        self.logger = logging.getLogger(__name__)

        # queue name (useful in Scheduler object when swapping queues)
        self.queue_name = queue_name

        # list of ObservingPrograms
        self.observing_programs = queue_configuration.build_observing_programs()

        # defaults to handle time-windowed queues
        self.is_TOO = False
        self.validity_window = None

        # Hack for greedy queues
        self.requests_in_window = True

        if 'validity_window_mjd' in queue_configuration.config:
            window = queue_configuration.config['validity_window_mjd']
            if window is not None:
                assert(len(window) == 2)
                self.set_validity_window_mjd(window[0], window[1])
            else:
                self.validity_window = None
        else:
            self.validity_window = None

        # flag to check if assign_nightly_requests has been called tonight
        self.queue_night = None

        # block on which the queue parameters were calculated
        self.queue_slot = None

        # number allowed requests by subprogram tonight 
        # (dict of (program_id, subprogram_name))
        self.requests_allowed = {}

        # the queue itself
        self.queue = pd.DataFrame()

        # should we only consider fields from one program in a given
        # observing block?  
        # CURRENTLY NOT IMPLEMENTED.
        self.block_programs = False

        if rp is None:
            # initialize an empty RequestPool
            self.rp = RequestPool()
        else:
            self.rp = rp

        if fields is None:
            self.fields = Fields()
        else:
            self.fields = fields

        self.missed_obs_queue = None

    def is_valid(self, time):
        if self.validity_window is None:
            return True

        window_start = self.validity_window[0]
        window_stop = self.validity_window[1]

        return window_start <= time <= window_stop

    def validity_window_mjd(self):
        if self.validity_window is None:
            return None

        return [self.validity_window[0].mjd, self.validity_window[1].mjd]

    def set_validity_window_mjd(self, window_start, window_stop):
        """Set the time at which this queue can run.

        Parameters
        ----------
        window_start : `float` 
            Modified Julian Date start time
        window_stop : `float` 
            Modified Julian Date end time
        """

        if window_start >= window_stop:
            raise ValueError("validity window start time must be less than end time")
        # rough checks
        if window_start <= Time('2017-01-01').mjd:
            raise ValueError(f"MJD likely out of range: {window_start}")
        if window_stop >= Time('2030-01-01').mjd:
            raise ValueError(f"MJD likely out of range: {window_stop}")

        self.validity_window = [Time(window_start,format='mjd'),
            Time(window_stop,format='mjd')]

    def compute_block_use(self):
        """Returns a dictionary with the fraction of blocks used by the queue,
        assuming observing starts at the beginning of the validity window"""
        

        if self.validity_window is None:
            raise ValueError('All blocks are valid')

        start_block = block_index(self.validity_window[0])
        obs_start_time = Time(self.validity_window[0],format='mjd')

        # greedy queues have no len until they have assignments made, so 
        # just use the validity window
        if len(self.queue) == 0:
            stop_block = block_index(self.validity_window[1])
            obs_end_time = self.validity_window[1]

        else:
            # with no weather, we start at the start of the window 
            if 'n_repeats' in self.queue.columns:
                n_obs = np.sum(self.queue.n_repeats)
                exp_time = np.sum(self.queue.exposure_time * self.queue.n_repeats)
            else:
                n_obs = len(self.queue)
                exp_time = np.sum(self.queue.exposure_time)
            obs_time = (exp_time * u.second) + n_obs * READOUT_TIME
            obs_end_time = self.validity_window[0] + obs_time

            stop_block = block_index(obs_end_time)
            # below breaks if the window is longer than the observations
            #stop_block = block_index(self.validity_window[1])

        assert obs_end_time > obs_start_time

        # compute fraction of the blocks used by the queue
        block_use = defaultdict(float)

        for block in np.arange(start_block, stop_block+1):

            block_use[block] = block_use_fraction(block, obs_start_time,
                                                  obs_end_time)

        return block_use

    def add_observing_program(self, observing_program):
        self.observing_programs.append(observing_program)

    def assign_nightly_requests(self, current_state, obs_log, 
            time_limit = 30 * u.second, block_use = defaultdict(float),
            timed_obs_count = defaultdict(int), skymaps = None):

        # clear previous request pool
        # missed obs and skymap greedy queues don't have observing programs
        if len(self.observing_programs) != 0:
            self.rp.clear_all_request_sets()
            # set number of allowed requests by program.
            self.determine_allowed_requests(current_state['current_time'],
                    obs_log, timed_obs_count = timed_obs_count)

        # can be used by field_selection_functions downstream
        program_fields = {}
        for program in self.observing_programs:
            key = (program.program_id, program.subprogram_name)
            program_fields[key] = \
                {'field_ids': program.field_ids,
                 'field_selection_function': program.field_selection_function,
                 'requests_allowed': self.requests_allowed[key]}

        for program in self.observing_programs:

            request_sets = program.assign_nightly_requests(
                current_state['current_time'], self.fields,
                obs_log, program_fields, block_programs=self.block_programs,
                skymaps = skymaps)
            for rs in request_sets:
                self.rp.add_request_sets(rs['program_id'], 
                            rs['subprogram_name'], rs['program_pi'],
                            rs['field_ids'], rs['filter_ids'], 
                            rs['intranight_gap'],
                            rs['exposure_time'],
                            rs['total_requests_tonight'])

#        assert(len(self.rp.pool) > 0)

        # any specific tasks needed)
        self._assign_nightly_requests(current_state, 
                time_limit = time_limit, block_use = block_use)

        # mark that we've set up the pool for tonight
        self.queue_night = np.floor(current_state['current_time'].mjd) 


    def adjust_program_exposures_tonight(self, obs_log, mjd_start, mjd_stop):
        """Use past history to adjust the number of exposures per program tonight.
        
        Counts exposures from the start of the month and equalizes any excess
        over NIGHTS_TO_REDISTRIBUTE or the number of nights to the end of 
        the month, whichever is less."""
        
        obs_count_by_program = obs_log.count_equivalent_obs_by_program(
                mjd_range = [mjd_start, mjd_stop])
        # drop engineering/commissioning
        obs_count_by_program = obs_count_by_program[
                obs_count_by_program['program_id'] != 0]
        obs_count_by_program.set_index('program_id', inplace=True)

        # if there are no observations, add zeros
        for program_id in PROGRAM_IDS:
            if program_id != 0:
                if program_id not in obs_count_by_program.index:
                    obs_count_by_program.loc[program_id] = 0

        total_obs = np.sum(obs_count_by_program['n_obs'])

        # infer the program fractions from the subprograms
        target_program_fractions = {propid:0 for propid in PROGRAM_IDS 
                if propid != 0}
        for op in self.observing_programs:
            target_program_fractions[op.program_id] = \
                    op.program_observing_time_fraction

        target_program_fractions = pd.Series(target_program_fractions) 
        target_program_fractions.index.name = 'program_id'
        target_program_fractions.name = 'target_fraction'

        target_program_nobs = target_program_fractions * total_obs
        target_program_nobs.name = 'target_program_nobs'

        # note that this gives 0 in case of no observations, as desired
        # have to do the subtraction backwords because of Series/DataFrame 
        # API nonsense
        delta_program_nobs = \
                -1*obs_count_by_program.subtract(target_program_nobs,
                    axis=0)

        NIGHTS_TO_REDISTRIBUTE = 5
        time = Time(mjd_stop,format='mjd')
        dtnow = time.to_datetime()
        if dtnow.month != 12:
            next_month_start_mjd = Time(datetime(dtnow.year,dtnow.month+1,1),
                    scale='utc').mjd
        else:
            next_month_start_mjd = Time(datetime(dtnow.year+1,1,1),
                    scale='utc').mjd
        nights_left_this_month = np.round(next_month_start_mjd - time.mjd)

        if nights_left_this_month > NIGHTS_TO_REDISTRIBUTE:
            divisor = NIGHTS_TO_REDISTRIBUTE
        else:
            divisor = nights_left_this_month
            if divisor == 0:
                divisor = 1

        delta_program_nobs /= divisor

        delta_program_nobs = np.round(delta_program_nobs).astype(int)

        return delta_program_nobs
        
    def adjust_subprogram_exposures_tonight(self, obs_log, mjd_start, mjd_stop):
        """Use past history to adjust the number of exposures per subprogram tonight.
        
        Counts exposures from the start of the month and equalizes any excess
        over NIGHTS_TO_REDISTRIBUTE or the number of nights to the end of 
        the month, whichever is less."""
        
        obs_count_by_subprogram_all = obs_log.count_equivalent_obs_by_subprogram(
                mjd_range = [mjd_start, mjd_stop])
        # drop engineering/commissioning
        obs_count_by_subprogram_all = obs_count_by_subprogram_all[
                obs_count_by_subprogram_all['program_id'] != 0]
        obs_count_by_subprogram_all.set_index(['program_id','subprogram_name'], 
                inplace=True)

        # only count the subprograms that are currently active.  This is
        # going to cause problems when the programs change--but we are going to
        # only use the subprogram balance for i-band
        obs_count_by_current_subprogram_dict = {}

        # if there are no observations, add zeros
        for op in self.observing_programs:
            idx = (op.program_id, op.subprogram_name) 
            if idx not in obs_count_by_subprogram_all.index:
                obs_count_by_current_subprogram_dict[idx] = 0
            else:
                obs_count_by_current_subprogram_dict[idx] = obs_count_by_subprogram_all.loc[idx,'n_obs']

        obs_count_by_subprogram = pd.Series(obs_count_by_current_subprogram_dict)
        obs_count_by_subprogram.name = 'n_obs'
        obs_count_by_subprogram.rename_axis(
                index=['program_id','subprogram_name'], inplace=True)

        total_obs = obs_count_by_subprogram.sum()

        # record the subprogram fractions
        target_subprogram_fractions = defaultdict(float)
        for op in self.observing_programs:
            target_subprogram_fractions[(op.program_id, op.subprogram_name)] = \
                    op.program_observing_time_fraction * op.subprogram_fraction

        target_subprogram_fractions = pd.Series(target_subprogram_fractions) 
#        target_program_fractions.index.name = 'program_id'
        target_subprogram_fractions.name = 'target_fraction'

        target_subprogram_nobs = target_subprogram_fractions * total_obs
        target_subprogram_nobs.name = 'target_subprogram_nobs'
        target_subprogram_nobs.rename_axis(
                index=['program_id','subprogram_name'], inplace=True)

        # note that this gives 0 in case of no observations, as desired
        # have to do the subtraction backwords because of Series/DataFrame 
        # API nonsense
        delta_subprogram_nobs = \
                -1*obs_count_by_subprogram.subtract(target_subprogram_nobs,
                    axis=0).fillna(0)

        NIGHTS_TO_REDISTRIBUTE = 5
        time = Time(mjd_stop,format='mjd')
        dtnow = time.to_datetime()
        if dtnow.month != 12:
            next_month_start_mjd = Time(datetime(dtnow.year,dtnow.month+1,1),
                    scale='utc').mjd
        else:
            next_month_start_mjd = Time(datetime(dtnow.year+1,1,1),
                    scale='utc').mjd
        nights_left_this_month = np.round(next_month_start_mjd - time.mjd)

        if nights_left_this_month > NIGHTS_TO_REDISTRIBUTE:
            divisor = NIGHTS_TO_REDISTRIBUTE
        else:
            divisor = nights_left_this_month
            if divisor == 0:
                divisor = 1

        delta_subprogram_nobs /= divisor

        delta_subprogram_nobs = np.round(delta_subprogram_nobs).astype(int)

        return delta_subprogram_nobs
        


    def determine_allowed_requests(self, time, obs_log, 
            timed_obs_count = defaultdict(int)):
        """Use count of past observations and expected observing time fractions
        to determine number of allowed requests tonight.
        
        Exclude observations already planned in timed queues."""

        self.requests_allowed = {}

        # rather than using equivalent obs, might be easier to work in 
        # exposure time directly?
        
        # enforce program balance on a monthly basis
        dtnow = time.to_datetime()
        month_start_mjd = Time(datetime(dtnow.year,dtnow.month,1),
                scale='utc').mjd
        
        delta_program_exposures_tonight = self.adjust_program_exposures_tonight(
            obs_log, month_start_mjd, time.mjd)
        # use this for i-band only
        delta_subprogram_exposures_tonight = self.adjust_subprogram_exposures_tonight(
            obs_log, month_start_mjd, time.mjd)
        
        self.logger.info(f'Change in allowed exposures: {delta_program_exposures_tonight}')
        self.logger.info(f'Needed change in allowed exposures by subprogram: {delta_subprogram_exposures_tonight}')
        self.logger.debug(f"Sum of change in allowed exposures by subprogram: {delta_subprogram_exposures_tonight.reset_index().groupby('program_id').agg(np.sum)}")
        self.logger.info(f'Number of timed observations: {timed_obs_count}')

        dark_time = approx_hours_of_darkness(time)
        
        # calculate subprogram fractions excluding list queues and TOOs
        scheduled_subprogram_sum = defaultdict(float)
        for op in self.observing_programs:
            # list queues and TOOs should set field_ids = [], but not None
            # OPs scheduled using field_selection_function will have 
            # field_ids = None
            if op.field_ids is not None:
                if len(op.field_ids) == 0:
                    continue
            scheduled_subprogram_sum[op.program_id] += \
                    op.subprogram_fraction

        for op in self.observing_programs:

            
            program_time_tonight = (
                dark_time * op.program_observing_time_fraction +  
                (delta_program_exposures_tonight.loc[op.program_id,'n_obs'] 
                - timed_obs_count[op.program_id]) * (EXPOSURE_TIME+READOUT_TIME))

            subprogram_time_tonight = (
                program_time_tonight * op.subprogram_fraction / 
                scheduled_subprogram_sum[op.program_id])

            n_requests = (subprogram_time_tonight.to(u.min) / 
                    op.time_per_exposure().to(u.min)).value[0]
            n_requests = np.round(n_requests).astype(int)

            # i_band program balance needs individual tuning due to 
            # longer cadence and filter blocking
            if op.subprogram_name == 'i_band':
                delta_i_nexp = delta_subprogram_exposures_tonight.loc[(2,'i_band')]
                if delta_i_nexp > 0:
                    self.logger.info(f'Adding {delta_i_nexp} additional i-band exposures')
                    n_requests += delta_i_nexp
                else:
                    self.logger.info(f'Implied change in i-band exposures is negative, skipping supplementation: {delta_i_nexp}')

            self.requests_allowed[(op.program_id, 
                op.subprogram_name)] = n_requests

        for key, n_requests in self.requests_allowed.items():
            if n_requests < 0:
                self.requests_allowed[key] = 0

        self.logger.info(self.requests_allowed)

    def next_obs(self, current_state, obs_log):
        """Given current state, return the parameters for the next request"""
        # don't store the telescope state locally!

        # check that assign_nightly_requests has been called tonight.
        if self.queue_type != 'list':
            if np.floor(current_state['current_time'].mjd) != self.queue_night:
                self.assign_nightly_requests(current_state, obs_log)

        # define functions that actually do the work in subclasses
        next_obs = self._next_obs(current_state, obs_log)

        # check if we have a disallowed observation, and reject it:
        if next_obs['target_limiting_mag'] < 0:
            self.logger.warning(f'Target is unobservable!  Removing from queue {next_obs}')
            self.remove_requests(next_obs['request_id'])
            next_obs = self.next_obs(current_state, obs_log)

        next_obs['queue_name'] = self.queue_name

        return next_obs

    def update_queue(self, current_state, obs_log, **kwargs):
        """Recalculate queue"""

        # define functions that actually do the work in subclasses
        return self._update_queue(current_state, obs_log)

    def remove_requests(self, request_id):
        """Remove a request from both the queue and the request set pool"""

        # define functions that actually do the work in subclasses
        return self._remove_requests(request_id)

    def move_program_to_missed_obs(self, program_id):
        """Remove all requests from the specified program_id to the missed obs queue."""

        # define functions that actually do the work in subclasses
        return self._move_program_to_missed_obs(program_id)

    def return_queue(self):
        """Return queue values, ordered in the expected sequence if possible"""

        queue = self._return_queue()

        cols = ['field_id','filter_id','exposure_time','program_id',
                'program_pi', 'max_airmass',
                'subprogram_name','ra','dec','ordered']
        if self.queue_type == 'gurobi':
            cols.append('slot_start_time')
        if self.queue_type == 'list':
            cols.append('mode_num')
            cols.append('ewr_num_images')
            cols.append('n_repeats')

        # pandas has gotten more picky if the columns aren't available
        out_cols = list(set(cols) & set(queue.columns))

        return queue.loc[:,out_cols]

class GurobiQueueManager(QueueManager):

    def __init__(self, queue_name, queue_configuration, **kwargs):
        super().__init__(queue_name, queue_configuration, **kwargs)
        self.block_obs_number = 0
        self.queue_type = 'gurobi'

    def _assign_nightly_requests(self, current_state, 
            time_limit = 30.*u.second, block_use = defaultdict(float)): 
        self._assign_slots(current_state, time_limit = time_limit, 
                block_use = block_use)

    def _next_obs(self, current_state, obs_log):
        """Select the highest value request."""

        # do the slot assignment at the beginning of the night 
        # (or if the queue is empty, which should be unusual)

        # if we've entered a new block, solve the TSP to sequence the requests
        if (block_index(current_state['current_time'])[0] != self.queue_slot):
            try:
                self._move_requests_to_missed_obs(self.queue_slot)
            except Exception as e:
                self.logger.exception(e)
                self.logger.error('Failed moving requests to missed obs!')
            self._sequence_requests_in_block(current_state)

        if (len(self.queue_order) == 0):
            raise QueueEmptyError("Ran out of observations this block.") 
        
        idx = self.queue_order[0]
        row = self.queue.loc[idx]
        if self.queue_slot in self.filter_by_slot:
            filter_id = int(self.filter_by_slot[self.queue_slot])
        else:
            raise QueueEmptyError("No requests in this slot!")

        next_obs = {'target_field_id': int(row['field_id']),
            'target_ra': row['ra'],
            'target_dec': row['dec'],
            'target_filter_id': filter_id,
            'target_program_id': int(row['program_id']),
            'target_subprogram_name': row['subprogram_name'],
            'target_program_pi': row['program_pi'],
            'target_exposure_time': row['exposure_time'] * u.second,
            'target_sky_brightness': 
                    self.block_sky_brightness.loc[idx,self.queue_slot][filter_id],
            'target_limiting_mag': 
                    self.block_lim_mags.loc[idx,self.queue_slot][filter_id],
            'target_metric_value':  
                    self.block_slot_metric.loc[idx,self.queue_slot][filter_id],
            'target_total_requests_tonight': int(row['total_requests_tonight']),
            'target_mode_num': 0,
            'target_num_images': 1,
            'request_id': idx}

#            'target_sky_brightness': self.queue.ix[idx].sky_brightness,
#            'target_limiting_mag': self.queue.ix[idx].limiting_mag,
#            'target_metric_value':  self.queue.ix[idx].value,
#            'target_request_number_tonight':

        return next_obs

    def _slot_metric(self, limiting_mag, dec):
        """Calculate metric for assigning fields to slots.

        penalizes volume for both extinction (airmass) and fwhm penalty
        due to atmospheric refraction, plus sky brightness from
        moon phase and distance
        == 1 for 21st mag.
        
        normalize metrics by maximum value at transit
        so low-declination fields are not penalized
        """
        #see 200430 notes

        metric = (10.**(0.6 * (limiting_mag - 21)) /
                (1-1e-4*(maximum_altitude(dec) - 90)**2.))
        # lock out -99 limiting mags even more aggressively
        return metric.where(limiting_mag > 0, -0.99)

    def _assign_slots(self, current_state, time_limit = 30*u.second, 
            block_use = defaultdict(float)):
        """Assign requests in the Pool to slots"""

        # check that the pool has fields in it
        if len(self.rp.pool) == 0:
            raise QueueEmptyError("No fields in pool")

        # join with fields so we have the information we need
        # make a copy so rp.pool and self.queue are not linked
        df = self.rp.pool.join(self.fields.fields, on='field_id').copy()

        # calculate limiting mag by block.  uses the block midpoint time
        blocks, times = nightly_blocks(current_state['current_time'], 
            time_block_size=TIME_BLOCK_SIZE)

        # remove the excluded blocks, if any.  Could do this in optimize.py
        # but it makes the optimization problem unneccesarily bigger
        # don't demand 100% of the block is used: tiny fractions lead to
        # infeasible models
        exclude_blocks = [b for (b,v) in block_use.items() if v > 0.95]

        self.logger.debug(f'Excluding completely filled blocks {exclude_blocks}')

        if len(exclude_blocks):
            cut_blocks = np.setdiff1d(blocks, exclude_blocks)
            cut_times = block_index_to_time(cut_blocks, 
                    current_state['current_time'], where='mid')
            blocks, times = cut_blocks, cut_times

        lim_mags = {}
        sky_brightnesses = {}
        decs = {}
        for bi, ti in zip(blocks, times):
            if 'altitude' in df.columns:
                df.drop('altitude', axis=1, inplace=True)
            if 'azimuth' in df.columns:
                df.drop('azimuth', axis=1, inplace=True)
            # use pre-computed blocks
            df_alt = self.fields.block_alt[bi]
            df_alt.name = 'altitude'
            df = df.join(df_alt, on='field_id')
            df_az = self.fields.block_az[bi]
            df_az.name = 'azimuth'
            df = df.join(df_az, on='field_id')
            for fid in FILTER_IDS:
                df_limmag, df_sky = \
                    compute_limiting_mag(df, ti, self.fields.Sky,
                            filter_id = fid)
                lim_mags[(bi, fid)] = df_limmag
                sky_brightnesses[(bi, fid)] = df_sky
                decs[(bi, fid)] = df.dec

        # this results in a MultiIndex on the *columns*: level 0 is block,
        # level 1 is filter_id.  df_metric.unstack() flattens it
        self.block_lim_mags = pd.DataFrame(lim_mags)
        self.block_sky_brightness = pd.DataFrame(sky_brightnesses)
        block_decs = pd.DataFrame(decs)
        self.block_slot_metric = self._slot_metric(self.block_lim_mags, 
                block_decs)

        # count the number of observations requested by filter
        df['n_reqs_tot'] = 0
        for fid in FILTER_IDS:
            df['n_reqs_{}'.format(fid)] = \
                df.filter_ids.apply(lambda x: np.sum([xi == fid for xi in x]))
            df['n_reqs_tot'] += df['n_reqs_{}'.format(fid)] 

        # prepare the data for input to gurobi
        #import shelve
        #s = shelve.open('tmp_vars.shelf')
        #s['block_lim_mags'] = self.block_lim_mags
        #s['block_slot_metric'] = self.block_slot_metric
        #s['df'] = df
        #s.close()

        self.request_sets_tonight, df_slots, dft = night_optimize(
            self.block_slot_metric, df, self.requests_allowed,
            time_limit = time_limit, block_use = block_use)

        grp = df_slots.groupby('slot')

        self.queued_requests_by_slot = grp['request_id'].apply(list)
        self.filter_by_slot = \
            grp['metric_filter_id'].apply(lambda x: np.unique(x)[0])

        # rework to dump output
        df_slots['scheduled'] = True
        dft.set_index(['request_id','slot','metric_filter_id'],inplace=True)
        df_slots.set_index(['request_id','slot','metric_filter_id'],inplace=True)
        dft = dft.join(df_slots,how='outer')
        dft['scheduled'] = dft['scheduled'].fillna(False)
        dft.reset_index(inplace=True)

        dft = pd.merge(dft,df[['field_id']],
            left_on='request_id', right_index=True)

        n_requests_scheduled = np.sum(dft['scheduled'])
        total_metric_value = np.sum(dft['scheduled']*dft['metric'])
        avg_metric_value = total_metric_value / n_requests_scheduled

        tot_avail_requests_bysubprogram = \
                df.groupby(['program_id','subprogram_name'])['n_reqs_tot'].agg(np.sum)
        tot_avail_requests_bysubprogram.name = 'available'

        # use self.requests_allowed and join this all up

        nscheduled_requests_bysubprogram = \
                dft.loc[dft['scheduled'],['program_id','subprogram_name']].groupby(['program_id','subprogram_name']).agg(len)
        nscheduled_requests_bysubprogram.name = 'scheduled'

        # reformat requests_allowed for joining
        mux = pd.MultiIndex.from_tuples(self.requests_allowed.keys(),
                names = ['program_id','subprogram_name'])
        df_allowed = pd.DataFrame(list(self.requests_allowed.values()),
                index=mux,columns=['allowed'])

        df_summary = df_allowed.join(tot_avail_requests_bysubprogram).join(nscheduled_requests_bysubprogram)
        self.logger.info(df_summary)

        self.logger.info(f'{n_requests_scheduled} requests scheduled')
        self.logger.info(f'{total_metric_value:.2f} total metric value; ' 
               f'{avg_metric_value:.2f} average per request')

        # this is not ideal for 
        tnow = current_state['current_time']
        yymmdd = tnow.iso.split()[0][2:].replace('-','')
        solution_outfile = f'{BASE_DIR}/../sims/gurobi_solution_{yymmdd}.csv'

        before_noon_utc = (tnow.mjd - np.floor(tnow.mjd)) < 0.5
        
        # avoid clobbering the solution file with restarts after observing has
        # completed
        if before_noon_utc or (not os.path.exists(solution_outfile)):
            dft.drop(columns=['Yrtf']).to_csv(solution_outfile)

    def _sequence_requests_in_block(self, current_state):
        """Solve the TSP for requests in this slot"""

        self.queue_slot = block_index(current_state['current_time'])[0]

        # raise an error if there are missing blocks--potentially due to
        # excluded blocks
        if self.queue_slot not in self.queued_requests_by_slot.index:
            raise QueueEmptyError(f"Current block {self.queue_slot} is not stored")

        # retrieve requests to be observed in this block
        req_list = self.queued_requests_by_slot.loc[self.queue_slot]

        # request_set ids should be unique per block
        assert( (len(set(req_list)) == len(req_list) ) )

        if np.all(np.isnan(req_list)):
            raise QueueEmptyError("No requests assigned to this block")

        idx = pd.Index(req_list)

        # reconstruct
        df = self.rp.pool.loc[idx].join(self.fields.fields, on='field_id').copy()
        az = self.fields.block_az[self.queue_slot]
        df = df.join(az, on='field_id')

        # now prepend the CALSTOW positoin so we can minimize slew from
        # filter exchanges 
        # Need to use current HA=0
        df_blockstart = pd.DataFrame({'ra':HA_to_RA(0,
            current_state['current_time']).to(u.degree).value,
            'dec':-48.,'azimuth':180.},index=[0])
        df_fakestart = pd.concat([df_blockstart,df],sort=True)

        # compute overhead time between all request pairs
        
        # compute pairwise slew times by axis for all pointings
        slews_by_axis = {}
        def coord_to_slewtime(coord, axis=None):
            c1, c2 = np.meshgrid(coord, coord)
            dangle = np.abs(c1 - c2)
            angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
            return slew_time(axis, angle * u.deg)

        slews_by_axis['dome'] = coord_to_slewtime(
            df_fakestart['azimuth'], axis='dome')
        slews_by_axis['dec'] = coord_to_slewtime(
            df_fakestart['dec'], axis='dec')
        slews_by_axis['ra'] = coord_to_slewtime(
            df_fakestart['ra'], axis='ha')

        maxradec = np.maximum(slews_by_axis['ra'], slews_by_axis['dec'])
        maxslews = np.maximum(slews_by_axis['dome'], maxradec)
        # impose a penalty on zero-length slews (which by construction
        # in this mode are from different programs)
        wnoslew = maxslews == 0
        maxslews[wnoslew] = READOUT_TIME * 10.
        overhead_time = np.maximum(maxslews, READOUT_TIME)

        tsp_order, tsp_overhead_time = tsp_optimize(overhead_time.value)

        # remove the fake starting point.  tsp_optimize always starts with
        # the first observation in df, which by construction is our fake point,
        # so we can simply cut it off.
        tsp_order = tsp_order[1:]
        assert(0 not in tsp_order)

        # tsp_order is 0-indexed from overhead time, so I need to
        # reconstruct the request_id
        self.queue_order = df_fakestart.index.values[tsp_order]
        self.queue = df

    def _move_requests_to_missed_obs(self, queue_slot, program_id=None):
        """After a block is expired, move any un-observed requests into the missed_obs queue."""
        #self.queue should have any remaining obs
        if len(self.queue):
            cols = ['program_id', 'subprogram_name', 'program_pi', 'field_id', 
                    'intranight_gap_min', 'exposure_time', 'probability']
            # it's a little confusing, because each queue entry has all of the
            # filter_ids from the original request set.  So we have to 
            # make a pool that only has single filters in it.
            filter_id = int(self.filter_by_slot[queue_slot])
            if program_id is not None:
                wp = self.queue['program_id'] == program_id
                if np.sum(wp):
                    missed_obs = self.queue.loc[wp,cols].copy()
                    self.logger.info(f"Moving {len(missed_obs)} requests (program {program_id}, filter {filter_id}) to the missed_obs queue: {missed_obs.loc[:,['subprogram_name','field_id']]}")
                else:
                    self.logger.info(f"No requests for program {program_id}, filter {filter_id} in current queue_slot {queue_slot}.")
                    return
            else:
                missed_obs = self.queue.loc[:,cols].copy()
            missed_obs['filter_ids'] = pd.Series([[filter_id] for i in missed_obs.index],index=missed_obs.index)
            missed_obs['total_requests_tonight'] = 1

            self.logger.info(f"Saving {len(missed_obs)} requests (filter {filter_id}) to the missed_obs queue: {missed_obs.loc[:,['subprogram_name','field_id']]}")

            # the missed obs RequestPool wants request *sets*, so find out
            # if previous requests were missed
            rows_to_append = []
            for idx, row in missed_obs.iterrows():
                if idx in self.missed_obs_queue.rp.pool.index:
                    assert(len(self.missed_obs_queue.rp.pool.loc[idx] == 1))
                    self.missed_obs_queue.rp.pool.loc[idx,'filter_ids'].append(filter_id) 
                    self.missed_obs_queue.rp.pool.loc[idx,'total_requests_tonight'] += 1
                else:
                    rows_to_append.append(row)
            self.missed_obs_queue.rp.pool = pd.concat([self.missed_obs_queue.rp.pool, pd.DataFrame(rows_to_append)])

        else:
            self.logger.debug(f'No remaining queued observations in slot {queue_slot}')

    def _remove_requests(self, request_set_id):
        """Remove a request from both the queue and the pool.
        
        Note that gurobi queue uses request_set_id to index."""

        # should be the topmost item
        assert (self.queue_order[0] == request_set_id)
        self.queue_order = self.queue_order[1:]
        row = self.queue.loc[request_set_id]
        self.queue = self.queue.drop(request_set_id)
        # (past slot assignments are still in self.queued_requests_by_slot)
        # (we will only reuse the RequestPool if we do recomputes)
        self.rp.remove_request(request_set_id, 
                self.filter_by_slot.loc[self.queue_slot])

    def _move_program_to_missed_obs(self, program_id):
        """Move all requests from program_id to missed_obs"""

        # move requests from current queue
        if self.queue_slot is not None:
            self._move_requests_to_missed_obs(self.queue_slot, program_id = program_id)

        # move request sets from the pool
        wpid = self.rp.pool['program_id'] == program_id
        if np.sum(wpid):
            self.missed_obs_queue.rp.pool = pd.concat([self.missed_obs_queue.rp.pool, 
                                                       self.rp.pool.loc[wpid]])
            self.rp.pool = self.rp.pool.loc[~wpid,:]

            # we also need to delete them self.queued_requests_by_slot
            slots = self.queued_requests_by_slot.index.values
            ids_to_remove = [di for di,dv in wpid.items() if dv]
            for slot in slots:
                self.queued_requests_by_slot.loc[slot] = [
                        req for req in self.queued_requests_by_slot.loc[slot]
                        if req not in ids_to_remove]

        self.logger.info(f'Moved {np.sum(wpid)} request sets from program id {program_id} to missed_obs')
            

    def _return_queue(self):

        # start by setting up the current slot
        if len(self.queue) > 0:
            queue = self.queue.loc[self.queue_order].copy()
            queue.loc[:,'ordered'] = True
            queue.loc[:,'slot_start_time'] = block_index_to_time(
                    self.queue_slot, Time.now(), where='start')[0].iso
        else:
            # before the night starts, the queue is empty
            queue = self.queue.copy()

        # now loop over upcoming slots, ensuring they are sorted (should be)
        slots = self.queued_requests_by_slot.index.values
        slots = np.sort(slots)

        for slot in slots:
            if (self.queue_slot is not None):
                if slot <= self.queue_slot:
                    continue
            slot_requests = self.queued_requests_by_slot.loc[slot]
 
            if len(slot_requests):
                idx = pd.Index(slot_requests)
                # reconstruct
                df = self.rp.pool.loc[idx].join(self.fields.fields, on='field_id').copy()
                df.loc[:,'filter_id'] = self.filter_by_slot[slot]
                df.loc[:,'ordered'] = False
                df.loc[:,'slot_start_time'] = block_index_to_time(slot,
                        Time.now(), where='start').iso[0]
                queue = pd.concat([queue,df])
            

        return queue


class GreedyQueueManager(QueueManager):

    def __init__(self, queue_name, queue_configuration, **kwargs):
        super().__init__(queue_name, queue_configuration, **kwargs)
        self.time_of_last_filter_change = None
        self.min_time_before_filter_change = TIME_BLOCK_SIZE
        self.queue_type = 'greedy'

    def _assign_nightly_requests(self, current_state,
            time_limit = 30.*u.second, block_use = defaultdict(float)):
        # initialize the time of last filter change
        if self.time_of_last_filter_change is None:
            self.time_of_last_filter_change = current_state['current_time']

    def _next_obs(self, current_state, obs_log):
        """Select the highest value request."""

        # since this is a greedy queue, we update the queue after each obs
        # for speed, only do the whole recalculation if we're in a new slot
#        if ((block_index(current_state['current_time'])[0] != self.queue_slot)
#                or (len(self.queue) == 0)):
#            self._update_queue(current_state)
#        else:
#            # otherwise just recalculate the overhead times
#            _ = self._update_overhead(current_state)

        # to get the "on the fly" cadence windows to work I have to 
        # run the whole queue every time right now...
        self._update_queue(current_state, obs_log)

        # in case this wasn't initialized by assign_nightly_requests
        if self.time_of_last_filter_change is None:
            self.time_of_last_filter_change = current_state['current_time']

        # check if filter changes are allowed yet
        if ((current_state['current_time'] - self.time_of_last_filter_change)
                < self.min_time_before_filter_change):
            # only consider observations in the current filter
            queue = self.queue[self.queue['filter_id'] == current_state['current_filter_id']]
            # unless there are no more observations, in which case allow a
            # change
            if len(queue) == 0:
                queue = self.queue
        else:
            # allow filter changes if desired
            queue = self.queue

        # request_id of the highest value request
        max_idx = queue.value.idxmax()
        row = queue.loc[max_idx]

        next_obs = {'target_field_id': int(row['field_id']),
            'target_ra': row['ra'],
            'target_dec': row['dec'],
            'target_filter_id': row['filter_id'],
            'target_program_id': row['program_id'],
            'target_subprogram_name': row['subprogram_name'],
            'target_program_pi': row['program_pi'],
            'target_exposure_time': row['exposure_time'] * u.second,
            'target_sky_brightness': row['sky_brightness'],
            'target_limiting_mag': row['limiting_mag'],
            'target_metric_value':  row['value'],
            'target_total_requests_tonight': row['total_requests_tonight'],
            'target_mode_num': 0,
            'target_num_images': 1,
            'request_id': max_idx}

        return next_obs

    def _metric(self, df):
        """Calculate metric for prioritizing fields.

        Penalizes volume for both extinction (airmass) and fwhm penalty
        due to atmospheric refraction, plus sky brightness from
        moon phase and distance, overhead time
        == 1 for 21st mag, 15 sec overhead.
        Normalize by value at transit."""

        return 10.**(0.6 * (df['limiting_mag'] - 21)) / \
            (1-1e-4*(maximum_altitude(df['dec']) - 90)**2.) / \
            ((EXPOSURE_TIME.value + df['overhead_time']) /
             (EXPOSURE_TIME.value + 10.)) * df['probability']

    def _update_overhead(self, current_state, df=None):
        """recalculate overhead values without regenerating whole queue"""

        inplace = df is None

        if inplace:
            # no dataframe supplied, so replace existing self.queue on exit
            df = self.queue
            df.drop(['overhead_time', 'altitude', 'azimuth'], axis=1,
                    inplace=True)

        # compute readout/slew overhead times, plus current alt/az
        df_overhead, df_altaz = self.fields.overhead_time(current_state)

        # nb: df has index request_id, not field_id
        df = pd.merge(df, df_overhead, left_on='field_id', right_index=True)
        df = pd.merge(df, df_altaz, left_on='field_id', right_index=True)
        df.rename(columns={'alt': 'altitude', 'az': 'azimuth'}, inplace=True)

        # add overhead for filter changes
        w = df['filter_id'] != current_state['current_filter_id']
        if np.sum(w):
            df.loc[w, 'overhead_time'] += FILTER_CHANGE_TIME.to(u.second).value

        if inplace:
            df.loc[:, 'value'] = self._metric(df)
            self.queue = df

        return df

    def _update_queue(self, current_state, obs_log):
        """Calculate greedy weighting of requests in the Pool using current
        telescope state only"""

        # store block index for which these values were calculated
        self.queue_slot = block_index(current_state['current_time'])[0]

        # check that the pool has fields in it
        if len(self.rp.pool) == 0:
            raise QueueEmptyError("No fields in pool")

        # join with fields so we have the information we need
        # make a copy so rp.pool and self.queue are not linked
        df_rs = self.rp.pool.join(self.fields.fields, on='field_id').copy()

        # now expand the dataframe of request sets to a dataframe with one
        # row per obs.  
        requests = []   
        for request_set_id, row in df_rs.iterrows():
            rdict = row.to_dict()
            filter_ids = rdict.pop('filter_ids')
            for filter_id in filter_ids:
                ri = rdict.copy()
                ri['filter_id'] = filter_id
                ri['request_set_id'] = request_set_id
                requests.append(ri)
        df = pd.DataFrame(requests)
            
        df = self._update_overhead(current_state, df=df)

        # start with conservative altitude cut;
        # airmass weighting applied naturally below
        # also make a copy because otherwise it retains knowledge of
        # (discarded) previous reference and raises SettingWithCopyWarnings
        df = df.loc[df['altitude'] > 20, :].copy()

        if len(df) == 0:
            raise QueueEmptyError("No fields in queue above altitude cut")

        # if restricting to one program per block, drop other programs
        if self.block_programs:
            current_block_program = PROGRAM_BLOCK_SEQUENCE[
                self.queue_slot % LEN_BLOCK_SEQUENCE]
            df = df.loc[df['program_id'] == current_block_program, :]

        cadence_cuts = enough_gap_since_last_obs(df,
            current_state,obs_log)

        self.requests_in_window = np.sum(cadence_cuts) > 0
        if ~self.requests_in_window:
            self.logger.warning(calc_queue_stats(df, current_state,
                intro="No fields with observable cadence windows.  Queue in progress:"))
            raise QueueEmptyError("No fields with observable cadence windows")
        # also make a copy because otherwise it retains knowledge of
        # (discarded) previous reference and raises SettingWithCopyWarnings
        df = df.loc[cadence_cuts, :].copy()

        # compute airmasses by field_id
        # airmass = zenith_angle_to_airmass(90. - df_alt)
        # airmass.name = 'airmass'
        # df = pd.merge(df, pd.DataFrame(airmass),
        #              left_on='field_id', right_index=True)
        # airmass cut (or add airmass weighting to value below)
        # df = df[(df['airmass'] <= MAX_AIRMASS) & (df['airmass'] > 0)]

        df_limmag, df_sky = compute_limiting_mag(df,
                current_state['current_time'], self.fields.Sky)
        df.loc[:, 'limiting_mag'] = df_limmag
        df.loc[:, 'sky_brightness'] = df_sky

        #df_limmag.name = 'limiting_mag'
        #df = pd.merge(df, df_limmag, left_on='field_id', right_index=True)

        df.loc[:, 'value'] = self._metric(df)

        self.queue = df

    def _remove_requests(self, request_id):
        """Remove a request from both the queue and the request pool"""

        row = self.queue.loc[request_id]

        self.queue = self.queue.drop(request_id)
        self.rp.remove_request(row['request_set_id'], row['filter_id'])

    def _move_program_to_missed_obs(self, program_id):
        raise NotImplementedError

    def _return_queue(self):

        if 'value' in self.queue.columns:
            queue = self.queue.sort_values('value',ascending=False).copy()
        else:
            queue = self.queue.copy()

        # we have put these in value order but the sequence can change
        queue['ordered'] = False
        return queue



class ListQueueManager(QueueManager):
    """Simple Queue that returns observations in order."""

    def __init__(self, queue_name, queue_configuration, fields=None, **kwargs):
        self.queue_type = 'list'

        # queue name (useful in Scheduler object when swapping queues)
        self.queue_name = queue_name

        if fields is None:
            self.fields = Fields()
        else:
            self.fields = fields

        # the queue itself
        self.load_list_queue(queue_configuration.config['targets'])

        if 'validity_window_mjd' in queue_configuration.config:
            window = queue_configuration.config['validity_window_mjd']
            if window is not None:
                assert(len(window) == 2)
                assert(window[1] > window[0])
                self.validity_window = [Time(window[0],format='mjd'),
                    Time(window[1],format='mjd')]
            else:
                self.validity_window = None
        else:
            self.validity_window = None
            
        self.is_TOO = queue_configuration.config['targets'][0]['subprogram_name'].startswith('ToO')

    def _assign_nightly_requests(self, current_state,
            **kwargs):
        pass

    def _update_queue(self, current_state, obs_log):
        pass

    def load_list_queue(self, queue_dict_list, append=False):
        """Initialize an ordered queue.

        queue_dict_list is a list of dicts, one per observation"""
        
        df = pd.DataFrame(queue_dict_list)

        # check that major columns are included
        required_columns = ['field_id','program_id', 'subprogram_name',
                'filter_id', 'program_pi']

        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f'Missing required column {col}')


        # by default use field ids alone to specify pointings, 
        # but allow manual ra/dec if needed
        if ('ra' not in df.columns) and ('dec' not in df.columns):
            queue = df.join(self.fields.fields, on='field_id', how='inner').sort_index().copy()
        else:
            queue = df

        # if some of the field ids are bad, there will be missing rows
        if len(queue) != len(df):
            raise ValueError('One or more field ids are malformed: {}'.format(
                df.index.difference(self.fields.fields.index)))

        # add standard keywords if not present
        if 'exposure_time' not in queue.columns:
            queue['exposure_time'] = EXPOSURE_TIME.to(u.second).value
        if 'max_airmass' not in queue.columns:
            queue['max_airmass'] = MAX_AIRMASS
        if 'n_repeats' not in queue.columns:
            queue['n_repeats'] = 1
        if 'mode_num' not in queue.columns:
            queue['mode_num'] = 0
        if 'ewr_num_images' not in queue.columns:
            queue['num_images'] = 1
        else:
            queue['num_images'] = queue['ewr_num_images']

        if append:
            self.queue = pd.concat([self.queue,queue],ignore_index=True)
        else:
            self.queue = queue

    def _next_obs(self, current_state, obs_log):
        """Return the next observation in the time ordered queue unless it has expired."""

        
        if len(self.queue) == 0:
            raise QueueEmptyError("No more observations in queue!")
        
        # take the next observation in line
        idx = 0

        while True:
            if idx == len(self.queue):
                raise QueueEmptyError("No valid observations in queue!")
            ra = self.queue.iloc[idx].ra
            ha = RA_to_HA(ra * u.degree, current_state['current_time']
                    ).to(u.degree).wrap_at(180.*u.degree).value
            dec = self.queue.iloc[idx].dec
            sc = coord.SkyCoord(ra,dec, unit=u.deg)
            airmass = altitude_to_airmass(
                    skycoord_to_altaz(sc, 
                        current_state['current_time']).alt.to(u.deg).value)
            if airmass >= self.queue.iloc[idx].max_airmass:
                idx += 1
                continue
            # Reed limits |HA| to < 5.95 hours (most relevant for circumpolar
            # fields not hit by the airmass cut)
            if np.abs(ha) >= (5.95 * u.hourangle).to(u.degree).value:
                idx += 1
                continue
            # 1) HA < -17.6 deg && Dec < -22 deg is rejected for both track & stow because of interference with FFI.
            if (ha <= -17.6) & (dec <= -22):
                idx += 1
                continue
             # West of HA -17.6 deg, Dec < -45 deg is rejected for tracking because of the service platform in the south.
            if (ha >= -17.6) & (dec <= -45):
                idx += 1
                continue
             # fabs(HA) > 3 deg is rejected for Dec < -46 to protect the shutter "ears".
            if (np.abs(ha) >= 3.) & (dec <= -46):
                idx += 1
                continue
             # dec > 87.5 is rejected
            if (dec > 87.5):
                idx += 1
                continue

            break

        
        next_obs = {'target_field_id': int(self.queue.iloc[idx].field_id),
            'target_ra': self.queue.iloc[idx].ra,
            'target_dec': self.queue.iloc[idx].dec,
            'target_filter_id': self.queue.iloc[idx].filter_id,
            'target_program_id': int(self.queue.iloc[idx].program_id),
            'target_subprogram_name': self.queue.iloc[idx].subprogram_name,
            'target_program_pi': self.queue.iloc[idx].program_pi,
            'target_exposure_time': self.queue.iloc[idx].exposure_time * u.second,
            'target_sky_brightness': 0.,
            'target_limiting_mag': 0.,
            'target_metric_value':  0.,
            'target_total_requests_tonight': 1,  
            'target_mode_num': int(self.queue.iloc[idx].mode_num),
            'target_num_images': int(self.queue.iloc[idx].num_images),
            'request_id': self.queue.index[idx]}

        return next_obs

    def _remove_requests(self, request_id):
        """Remove a request from the queue"""

        try:
            if self.queue.loc[request_id,'n_repeats'] > 1:
                self.queue.loc[request_id,'n_repeats'] -= 1
            else:    
                self.queue = self.queue.drop(request_id)
        except Exception:
            self.logger.exception(f'Failure removing request {request_id}')

    def _move_program_to_missed_obs(self, program_id):
        raise NotImplementedError

    def _return_queue(self):

        # by construction the list queue is already in order
        queue = self.queue.copy()
        queue['ordered'] = True
        return queue

class RequestPool(object):

    def __init__(self):
        # initialize empty dataframe to add to
        self.pool = pd.DataFrame()
        pass

    def add_request_sets(self, program_id, subprogram_name, program_pi,
                field_ids, filter_ids, intranight_gap, exposure_time, 
                total_requests_tonight, probability=1):
        """program_ids must be scalar"""

        assert (scalar_len(program_id) == 1) 
        assert (scalar_len(subprogram_name) == 1) 

        n_fields = scalar_len(field_ids)
        if n_fields == 1:
            # see if it's iterable or not
            try:
                iterator = iter(field_ids)
            except TypeError:
                # if not, assume it's a scalar and wrap in a list
                field_ids = [field_ids]

        n_probs = scalar_len(probability)
        if n_probs == 1:
            # blow it up to match field_ids
            probabilities = np.ones(n_fields) * probability
        else:
            assert(n_probs == n_fields)
            probabilities = probability

        # build df as a list of dicts
        request_sets = []
        for i, (field_id, prob_i)  in enumerate(zip(field_ids, probabilities)):
            request_sets.append({
                'program_id': program_id,
                'subprogram_name': subprogram_name,
                'program_pi': program_pi,
                'field_id': field_id,
                'filter_ids': filter_ids.copy(),
                # pandas doesn't play well with astropy quantities, so change
                # back to seconds
                'intranight_gap_min': intranight_gap.to(u.minute).value,
                'exposure_time': exposure_time.to(u.second).value,
                'total_requests_tonight': total_requests_tonight,
                'probability': prob_i})

        self.pool = pd.concat([self.pool, pd.DataFrame(request_sets)], 
            ignore_index=True)

    def n_request_sets(self):
        return len(self.pool)

    def remove_request_sets(self, request_set_ids):
        """Remove completed or otherwise unwanted requests by request_id

        request_ids : scalar or list
            requests to drop (index of self.pool)"""
        self.pool = self.pool.drop(request_set_ids)

    def remove_request(self, request_set_id, filter_id):
        """Remove single completed request from a request set. 

        request_set_id: scalar 
            request set to modify (index of self.pool)
        filter_id: scalar
            filter_id of completed observation"""

        rs = self.pool.loc[request_set_id].copy()
        filters = rs['filter_ids']
        # this is another step that shouldn't be necessary...
        filters.remove(filter_id)
        if len(filters) == 0:
            self.remove_request_sets(request_set_id)
        else:
            self.pool.at[request_set_id, 'filter_ids'] =  filters

    def clear_all_request_sets(self):
        self.pool = pd.DataFrame()


# utils for examining inputs

def calc_pool_stats(df, intro=""):
    """

    df = Q.rp.pool"""

    stats_str = intro + "\n"
    stats_str += "\t{} request sets\n".format(len(df))
    stats_str += "\t{} unique fields\n".format(len(set(df.field_id)))
    for prog_id in PROGRAM_IDS:
        w = df.program_id == prog_id
        stats_str += "\tProgram {}:\n".format(prog_id)
        stats_str += "\t\t{} request sets\n".format(np.sum(w))
        stats_str += "\t\t{} unique fields\n".format(
            len(set(df.loc[w, 'field_id'])))
        stats_str += "\t\t{} median requests tonight per field\n".format(
            np.median(df.loc[w, 'total_requests_tonight']))

    return stats_str


def calc_queue_stats(df, current_state, intro=""):
    """

    df = Q.queue"""

    stats_str = intro + "\n"
    stats_str += "\t{} queued requests\n".format(len(df))
    stats_str += "\t{} unique fields\n".format(len(set(df.field_id)))
    for prog_id in PROGRAM_IDS:
        w = df.program_id == prog_id
        stats_str += "\tProgram {}:\n".format(prog_id)

        if np.sum(w) == 0:
            stats_str += "\t\tNo queued requests!\n"
            continue

        stats_str += "\t\t{} requests\n".format(np.sum(w))
        stats_str += "\t\t{} unique fields\n".format(
            len(set(df.loc[w, 'field_id'])))
        walt = w & (df.loc[w, 'altitude'] > 20)
        stats_str += "\t\t{} fields above altitude cut\n".format(
            np.sum(walt))
#        wfirst = walt & (df.loc[walt, 'request_number_tonight'] == 1)
#        stats_str += "\t\t{} requests awaiting first obs tonight\n".format(
#            np.sum(wfirst))

    return stats_str
