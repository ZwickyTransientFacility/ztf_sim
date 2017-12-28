from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import object
from collections import defaultdict
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time
import astroplan
from .Fields import Fields
from .SkyBrightness import SkyBrightness
from .magnitudes import limiting_mag
from .optimize import request_set_optimize, slot_optimize, tsp_optimize
from .cadence import enough_gap_since_last_obs
from .constants import P48_loc, PROGRAM_IDS, FILTER_IDS, TIME_BLOCK_SIZE
from .constants import EXPOSURE_TIME, READOUT_TIME, FILTER_CHANGE_TIME, slew_time
from .constants import PROGRAM_BLOCK_SEQUENCE, LEN_BLOCK_SEQUENCE, MAX_AIRMASS
from .utils import skycoord_to_altaz, seeing_at_pointing
from .utils import altitude_to_airmass, airmass_to_altitude, RA_to_HA
from .utils import scalar_len, nightly_blocks, block_index, block_index_to_time

class QueueEmptyError(Exception):
    """Error class for when the nightly queue has no more fields"""
    pass


class QueueManager(object):

    def __init__(self, observing_programs=[], rp=None, fields=None,
                 block_programs=False, queue_name = 'default'):

        # queue name (useful in Scheduler object when swapping queues)
        self.queue_name = queue_name

        # list of ObservingPrograms
        self.observing_programs = observing_programs

        # block on which the queue parameters were calculated
        self.queue_slot = None

        # number allowed requests by subprogram tonight 
        # (dict of (program_id, subprogram_name))
        self.requests_allowed = {}

        # the queue itself
        self.queue = pd.DataFrame()

        # should we only consider fields from one program in a given
        # observing block?
        self.block_programs = block_programs

        if rp is None:
            # initialize an empty RequestPool
            self.rp = RequestPool()
        else:
            self.rp = rp

        if fields is None:
            self.fields = Fields()
        else:
            self.fields = fields

        self.Sky = SkyBrightness()

    def add_observing_program(self, observing_program):
        self.observing_programs.append(observing_program)

    def assign_nightly_requests(self, current_state, obs_log):
        # clear previous request pool
        self.rp.clear_all_request_sets()
        # set number of allowed requests by program.
        self.determine_allowed_requests(current_state['current_time'],
                obs_log)

        for program in self.observing_programs:

            request_sets = program.assign_nightly_requests(
                current_state['current_time'], self.fields,
                block_programs=self.block_programs)
            for rs in request_sets:
                self.rp.add_request_sets(rs['program_id'], 
                            rs['subprogram_name'], rs['field_ids'],
                            rs['filter_ids'], rs['total_requests_tonight'])

        assert(len(self.rp.pool) > 0)

        # any specific tasks needed)
        self._assign_nightly_requests(current_state)

    def determine_allowed_requests(self, time, obs_log):
        """Use count of past observations and expected observing time fractions
        to determine number of allowed requests tonight."""

        self.requests_allowed = defaultdict(int)
        
        obs_count_by_subprogram = obs_log.count_total_obs_by_subprogram()
        total_obs = np.sum(list(obs_count_by_subprogram.values()))
        for program in self.observing_programs:
            n_requests = program.number_of_allowed_requests(time)
            delta = np.round(
                (obs_count_by_subprogram[(program.program_id, 
                    program.subprogram_name)] -
                 program.program_observing_time_fraction * total_obs)
                 * program.subprogram_fraction)

            # TODO: tweak as needed
            # how quickly do we want to take to reach equalization?
            CATCHUP_FACTOR = 0.20
            n_requests -= np.round(delta * CATCHUP_FACTOR).astype(np.int)
            self.requests_allowed[(program.program_id, 
                program.subprogram_name)] += n_requests

        for key, n_requests in self.requests_allowed.items():
            if n_requests < 0:
                self.requests_allowed[key] = 0

    def next_obs(self, current_state):
        """Given current state, return the parameters for the next request"""
        # don't store the telescope state locally!

        # define functions that actually do the work in subclasses
        return self._next_obs(current_state)

    def update_queue(self, current_state, **kwargs):
        """Recalculate queue"""

        # define functions that actually do the work in subclasses
        return self._update_queue(current_state)

    def remove_requests(self, request_id):
        """Remove a request from both the queue and the request set pool"""

        # define functions that actually do the work in subclasses
        return self._remove_requests(request_id)

    def compute_limiting_mag(self, df, time, filter_id=None):
        """compute limiting magnitude based on sky brightness and seeing"""

        # copy df so we can edit the filter id if desired
        if filter_id is not None:
            df = df.copy()
            df['filter_id'] = filter_id

        # compute inputs for sky brightness
        sc = coord.SkyCoord(df['ra'], df['dec'], frame='icrs', unit='deg')
        sun = coord.get_sun(time)
        sun_altaz = skycoord_to_altaz(sun, time)
        moon = coord.get_moon(time, location=P48_loc)
        moon_altaz = skycoord_to_altaz(moon, time)
        df.loc[:, 'moonillf'] = astroplan.moon.moon_illumination(time)
        
        # WORKING AROUND BUG in moon distance!!!!  171110
        df.loc[:, 'moon_dist'] = moon.separation(sc).to(u.deg).value
        df.loc[:, 'moonalt'] = moon_altaz.alt.to(u.deg).value
        df.loc[:, 'sunalt'] = sun_altaz.alt.to(u.deg).value

        # check if the sun is up anywhere and break things if it isn't
        if np.sum(df['sunalt'] > -6) != 0:
            raise ValueError('Some pointings outside six-degree twilight!')

        # compute sky brightness
        # only have values for reasonable altitudes (set by R20_absorbed...)
        wup = df['altitude'] >= airmass_to_altitude(MAX_AIRMASS) 
        df.loc[wup, 'sky_brightness'] = self.Sky.predict(df[wup])

        # compute seeing at each pointing
        # TODO: allow variable zenith seeing?  or by band?
        df.loc[wup, 'seeing'] = seeing_at_pointing(2.0*u.arcsec, 
            df.loc[wup,'altitude'])

        df.loc[wup, 'limiting_mag'] = limiting_mag(EXPOSURE_TIME, 
            df.loc[wup, 'seeing'],
            df.loc[wup, 'sky_brightness'],
            filter_id = df.loc[wup,'filter_id'],
            altitude = df.loc[wup,'altitude'], SNR=5.)

        # renormalize limiting mags to the R-band range so we maintain 
        # the causal structure with airmass, etc. but can get i-band scheduled
        
        # bright time limiting mags (from PTF-trained model--see 170930 notes
        # and plot_sky_brightness_model.ipynb)
        mlim_bright_g = 19.9
        mlim_bright_r = 20.1
        mlim_bright_i = 19.5
        dm_g = (21.9-19.9)
        dm_r = (21.5-20.1)
        dm_i = (20.9-19.5)

        wg = df['filter_id'] == 1
        if np.sum(wg):
            df.loc[wg,'limiting_mag'] = \
                (df.loc[wg,'limiting_mag'] - mlim_bright_g) * dm_r/dm_g \
                + mlim_bright_r

        wi = df['filter_id'] == 3
        if np.sum(wi):
            df.loc[wi,'limiting_mag'] = \
                (df.loc[wi,'limiting_mag'] - mlim_bright_i) * dm_r/dm_i \
                + mlim_bright_r

        # assign a very bright limiting mag to the fields that are down 
        # so the metric goes to zero
        df.loc[~wup, 'limiting_mag'] = -99

        # assign a very bright limiting mag to the fields within 20 degrees of
        # the moon 
        wmoon = df['moon_dist'] < 20
        df.loc[wmoon, 'limiting_mag'] = -99

        ha_vals = RA_to_HA(df['ra'].values*u.degree, time)
        # for limits below, need ha-180-180
        ha_vals = ha_vals.wrap_at(180.*u.degree)
        ha = pd.Series(ha_vals.to(u.degree), index=df.index, name='ha')

        # lock out TCS limits
        
        # Reed limits |HA| to < 5.95 hours (most relevant for circumpolar
        # fields not hit by the airmass cut)
        whalimit = np.abs(ha) >= (5.95 * u.hourangle).to(u.degree).value
        df.loc[whalimit, 'limiting_mag'] = -99
        
        # 1) HA < -17.6 deg && Dec < -22 deg is rejected for both track & stow because of interference with FFI.
        
        w1 = (ha <= -17.6) & (df['dec'] <= -22)
        df.loc[w1, 'limiting_mag'] = -99

        # West of HA -17.6 deg, Dec < -45 deg is rejected for tracking because of the service platform in the south.  
        w2 = (ha >= -17.6) & (df['dec'] <= -45)
        df.loc[w2, 'limiting_mag'] = -99

        # fabs(HA) > 3 deg is rejected for Dec < -46 to protect the shutter "ears".  
        w3 = (np.abs(ha) >= 3.) & (df['dec'] <= -46)
        df.loc[w3, 'limiting_mag'] = -99

        return df['limiting_mag'], df['sky_brightness']

    def return_queue(self):
        """Return queue values, ordered in the expected sequence if possible"""

        queue = self._return_queue()

        cols = ['field_id','filter_id','exposure_time','program_id',
                'subprogram_name','ra','dec','ordered']
        if self.queue_type == 'gurobi':
            cols.append('slot_start_time')

        return queue.loc[:,cols]



class GurobiQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GurobiQueueManager, self).__init__(**kwargs)
        self.block_obs_number = 0
        self.queue_type = 'gurobi'

    def _assign_nightly_requests(self, current_state):
        self._assign_slots(current_state)

    def _next_obs(self, current_state):
        """Select the highest value request."""

        # do the slot assignment at the beginning of the night 
        # (or if the queue is empty, which should be unusual)
        # TODO: consider how to detect and handle weather losses--recompute?
        #if (len(self.queue) == 0):
        #    self._assign_slots(current_state)

        # if we've entered a new block, solve the TSP to sequence the requests
        if (block_index(current_state['current_time'])[0] != self.queue_slot):
            self._sequence_requests_in_block(current_state)


        # TODO: should now be able to simply pick up the next observation
        
        if len(self.queue_order) == 0:
            raise QueueEmptyError("Ran out of observations this block.") 
        
        idx = self.queue_order[0]
        row = self.queue.loc[idx]
        filter_id = int(self.filter_by_slot[self.queue_slot])

        # TODO: make the queue have the right datatypes
        next_obs = {'target_field_id': int(row['field_id']),
            'target_ra': row['ra'],
            'target_dec': row['dec'],
            'target_filter_id': filter_id,
            'target_program_id': int(row['program_id']),
            'target_subprogram_name': row['subprogram_name'],
            'target_exposure_time': EXPOSURE_TIME,
            'target_sky_brightness': 
                    self.block_sky_brightness.loc[idx,self.queue_slot][filter_id],
            'target_limiting_mag': 
                    self.block_lim_mags.loc[idx,self.queue_slot][filter_id],
            'target_metric_value':  
                    self.block_slot_metric.loc[idx,self.queue_slot][filter_id],
            'target_total_requests_tonight': int(row['total_requests_tonight']),
            'request_id': idx}

#            'target_sky_brightness': self.queue.ix[idx].sky_brightness,
#            'target_limiting_mag': self.queue.ix[idx].limiting_mag,
#            'target_metric_value':  self.queue.ix[idx].value,
#            'target_request_number_tonight':

        return next_obs

    def _slot_metric(self, limiting_mag):
        """Calculate metric for assigning fields to slots.

        penalizes volume for both extinction (airmass) and fwhm penalty
        due to atmospheric refraction, plus sky brightness from
        moon phase and distance
        == 1 for 21st mag."""
        return 10.**(0.6 * (limiting_mag - 21)) 

    def _assign_slots(self, current_state):
        """Assign requests in the Pool to slots"""

        # check that the pool has fields in it
        if len(self.rp.pool) == 0:
            raise QueueEmptyError("No fields in pool")

        # join with fields so we have the information we need
        # make a copy so rp.pool and self.queue are not linked
        df = self.rp.pool.join(self.fields.fields, on='field_id').copy()

        # calculate limiting mag by block 
        # TODO: and by filter
        blocks, times = nightly_blocks(current_state['current_time'], 
            time_block_size=TIME_BLOCK_SIZE)

        lim_mags = {}
        sky_brightnesses = {}
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
            # TODO: only using r-band right now
            for fid in FILTER_IDS:
                df_limmag, df_sky = \
                    self.compute_limiting_mag(df, ti, filter_id = fid)
                lim_mags[(bi, fid)] = df_limmag
                sky_brightnesses[(bi, fid)] = df_sky

        # this results in a MultiIndex on the *columns*: level 0 is block,
        # level 1 is filter_id.  df_metric.unstack() flattens it
        self.block_lim_mags = pd.DataFrame(lim_mags)
        self.block_sky_brightness = pd.DataFrame(sky_brightnesses)
        self.block_slot_metric = self._slot_metric(self.block_lim_mags)

        # count the number of observations requested by filter
        for fid in FILTER_IDS:
            df['n_reqs_{}'.format(fid)] = \
                df.filter_ids.apply(lambda x: np.sum([xi == fid for xi in x]))

        # prepare the data for input to gurobi
        #import shelve
        #s = shelve.open('tmp_vars.shelf')
        #s['block_lim_mags'] = self.block_lim_mags
        #s['block_slot_metric'] = self.block_slot_metric
        #s['df'] = df
        #s.close()

        # select request_sets for the night
        self.request_sets_tonight = request_set_optimize(
            self.block_slot_metric, df, self.requests_allowed)

        # optimize assignment into slots
        df_slots = slot_optimize(
            self.block_slot_metric.loc[self.request_sets_tonight], 
            df.loc[self.request_sets_tonight], self.requests_allowed)

        grp = df_slots.groupby('slot')

        self.queued_requests_by_slot = grp['request_id'].apply(list)
        self.filter_by_slot = \
            grp['metric_filter_id'].apply(lambda x: np.unique(x)[0])



    def _sequence_requests_in_block(self, current_state):
        """Solve the TSP for requests in this slot"""

        self.queue_slot = block_index(current_state['current_time'])[0]

        assert(self.queue_slot in self.queued_requests_by_slot.index)

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

        # compute overhead time between all request pairs
        
        # compute pairwise slew times by axis for all pointings
        slews_by_axis = {}
        def coord_to_slewtime(coord, axis=None):
            c1, c2 = np.meshgrid(coord, coord)
            dangle = np.abs(c1 - c2)
            angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
            return slew_time(axis, angle * u.deg)

        slews_by_axis['dome'] = coord_to_slewtime(
            df['azimuth'], axis='dome')
        slews_by_axis['dec'] = coord_to_slewtime(
            df['dec'], axis='dec')
        slews_by_axis['ra'] = coord_to_slewtime(
            df['ra'], axis='ha')

        maxradec = np.maximum(slews_by_axis['ra'], slews_by_axis['dec'])
        maxslews = np.maximum(slews_by_axis['dome'], maxradec)
        overhead_time = np.maximum(maxslews, READOUT_TIME)

        tsp_order, tsp_overhead_time = tsp_optimize(overhead_time.value)

        # tsp_order is 0-indexed from overhead time, so I need to
        # reconstruct the request_id
        self.queue_order = df.index.values[tsp_order]
        self.queue = df

        scheduled_time = len(tsp_order) * EXPOSURE_TIME + \
            tsp_overhead_time*u.second

        # TODO: some sort of monitoring of expected vs. block time vs. actual

    def _remove_requests(self, request_set_id, completed=False):
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

    def _return_queue(self):

        # start by setting up the current slot
        if len(self.queue) > 0:
            queue = self.queue.loc[self.queue_order].copy()
            queue.loc[:,'ordered'] = True
            queue.loc[:,'slot_start_time'] = block_index_to_time(
                    self.queue_slot, Time.now(), where='start').iso
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
 
            idx = pd.Index(slot_requests)
            # reconstruct
            df = self.rp.pool.loc[idx].join(self.fields.fields, on='field_id').copy()
            df.loc[:,'filter_id'] = self.filter_by_slot[slot]
            df.loc[:,'ordered'] = False
            df.loc[:,'slot_start_time'] = block_index_to_time(slot,
                    Time.now(), where='start').iso
            queue = queue.append(df)
        

        return queue


class GreedyQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GreedyQueueManager, self).__init__(**kwargs)
        self.time_of_last_filter_change = None
        self.min_time_before_filter_change = TIME_BLOCK_SIZE
        self.queue_type = 'greedy'

    def _assign_nightly_requests(self, current_state):
        # initialize the time of last filter change
        if self.time_of_last_filter_change is None:
            self.time_of_last_filter_change = current_state['current_time']

    def _next_obs(self, current_state):
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
        self._update_queue(current_state)

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
        max_idx = queue.value.argmax()
        row = queue.loc[max_idx]


        # enforce program balance
        progid = row['program_id']
        subprogram = row['subprogram_name']
        # TODO: use obs_log to query for requests completed by this subprogram
        # tonight
        #if ((self.requests_completed[(progid, subprogram)] + 1) >= 
        #        self.requests_allowed[(progid, subprogram)]):
        if False:
            # this program has used up all its obs for tonight.
            # remove all requests for this program from the pool 
            w = ((self.rp.pool.program_id == progid) &
                 (self.rp.pool.subprogram_name == subprogram))
            self.rp.remove_requests(self.rp.pool[w].index.tolist())
            # reset the queue
            self._update_queue(current_state)
            # and request a new next_obs
            next_obs = self._next_obs(current_state)
        else:
            next_obs = {'target_field_id': row['field_id'],
                'target_ra': row['ra'],
                'target_dec': row['dec'],
                'target_filter_id': row['filter_id'],
                'target_program_id': row['program_id'],
                'target_subprogram_name': row['subprogram_name'],
                'target_exposure_time': EXPOSURE_TIME,
                'target_sky_brightness': row['sky_brightness'],
                'target_limiting_mag': row['limiting_mag'],
                'target_metric_value':  row['value'],
                'target_total_requests_tonight': row['total_requests_tonight'],
                'request_id': max_idx}

        return next_obs

    def _metric(self, df):
        """Calculate metric for prioritizing fields.

        penalizes volume for both extinction (airmass) and fwhm penalty
        due to atmospheric refraction, plus sky brightness from
        moon phase and distance, overhead time
        == 1 for 21st mag, 15 sec overhead."""
        # df.loc[:, 'value'] =
        return 10.**(0.6 * (df['limiting_mag'] - 21)) / \
            ((EXPOSURE_TIME.value + df['overhead_time']) /
             (EXPOSURE_TIME.value + 15.))

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
        # TODO: standardize this naming
        df.rename(columns={'alt': 'altitude', 'az': 'azimuth'}, inplace=True)

        # add overhead for filter changes
        w = df['filter_id'] != current_state['current_filter_id']
        if np.sum(w):
            df.loc[w, 'overhead_time'] += FILTER_CHANGE_TIME.to(u.second).value

        if inplace:
            df.loc[:, 'value'] = self._metric(df)
            self.queue = df

        return df

    def _update_queue(self, current_state):
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



        # use cadence functions to compute requests with active cadence windows
        # this is slow, so do it after our altitude cut
        # TODO: consider instead
        # df.apply(time_since_obs,args=(current_state,),axis=1)
        #in_window = {}
        #for idx, row in df.iterrows():
        #    # this is way, way slower
        #    # df['in_cadence_window'].ix[idx] = \
        #    in_window[idx] = eval(
        #        '{}(row, current_state)'.format(row['cadence_func']))
        #cadence_cuts = pd.Series(in_window)

        cadence_cuts  = df.apply(enough_gap_since_last_obs, 
            args=(current_state,),axis=1)
        

        # TODO: handle if cadence cuts returns no fields
        if np.sum(cadence_cuts) == 0:
            print(calc_queue_stats(df, current_state,
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

        df_limmag, df_sky = self.compute_limiting_mag(df,
                current_state['current_time'])
        df.loc[:, 'limiting_mag'] = df_limmag
        df.loc[:, 'sky_brightness'] = df_sky

        #df_limmag.name = 'limiting_mag'
        #df = pd.merge(df, df_limmag, left_on='field_id', right_index=True)

        df.loc[:, 'value'] = self._metric(df)

        self.queue = df

    def _remove_requests(self, request_id)
        """Remove a request from both the queue and the request pool"""

        row = self.queue.loc[request_id]

        self.queue = self.queue.drop(request_id)
        self.rp.remove_request(row['request_set_id'], row['filter_id'])

    def _return_queue(self):

        queue = self.queue.sort_values('value',ascending=False).copy()
        # we have put these in value order but the sequence can change
        queue['ordered'] = False
        return queue



class ListQueueManager(QueueManager):
    """Simple Queue that returns observations in order."""

    def __init__(self, **kwargs):
        super(ListQueueManager, self).__init__(**kwargs)
        self.queue_type = 'list'

    def _assign_nightly_requests(self, current_state):
        raise NotImplementedError("ListQueueManager should be loaded by load_queue()")

    def _update_queue(self, current_state):
        pass

    def load_list_queue(self, queue_dict_list, append=False):
        """Initialize an ordered queue.

        queue_dict_list is a list of dicts, one per observation"""
        
        df = pd.DataFrame(queue_dict_list)

        # check that major columns are included
        required_columns = ['field_id','program_id', 'subprogram_name',
                'filter_id']

        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f'Missing required column {col}')


        # by default use field ids alone to specify pointings, 
        # but allow manual ra/dec if needed
        if ('ra' not in df.columns) and ('dec' not in df.columns):
            queue = df.join(self.fields.fields, on='field_id', how='inner').copy()

        # if some of the field ids are bad, there will be missing rows
        if len(queue) != len(df):
            raise ValueError('One or more field ids are malformed: {}'.format(
                df.index.difference(self.fields.fields.index)))

        # add standard keywords if not present
        if 'exposure_time' not in queue.columns:
            queue['exposure_time'] = EXPOSURE_TIME
        if 'max_airmass' not in queue.columns:
            queue['max_airmass'] = MAX_AIRMASS
        if 'n_repeats' not in queue.columns:
            queue['n_repeats'] = 1


        if append:
            self.queue = self.queue.append(queue, ignore_index=True)
        else:
            self.queue = queue



    def _next_obs(self, current_state):
        """Return the next observation in the time ordered queue unless it has expired."""

        
        if len(self.queue) == 0:
            raise QueueEmptyError("No more observations in queue!")
        
        # take the next observation in line
        idx = 0

        # TODO: check if it's past the max airmass
        while True:
            if len(self.queue) == 0:
                raise QueueEmptyError("No more observations in queue!")
            ra = self.queue.iloc[idx].ra
            dec = self.queue.iloc[idx].dec
            sc = coord.SkyCoord(ra,dec, unit=u.deg)
            airmass = altitude_to_airmass(
                    skycoord_to_altaz(sc, 
                        current_state['current_time']).alt.to(u.deg).value)
            if airmass <= self.queue.iloc[idx].max_airmass:
                break
            else:
                print('Above max airmass, removing from queue!')
                self._remove_requests(self.queue.index[idx])
                # this effectively pops off iloc=0, so we can keep idx=0.  

        
        # TODO: think if we want a time-based expiration as well

        next_obs = {'target_field_id': int(self.queue.iloc[idx].field_id),
            'target_ra': self.queue.iloc[idx].ra,
            'target_dec': self.queue.iloc[idx].dec,
            'target_filter_id': self.queue.iloc[idx].filter_id,
            'target_program_id': int(self.queue.iloc[idx].program_id),
            'target_subprogram_name': self.queue.iloc[idx].subprogram_name,
            'target_exposure_time': self.queue.iloc[idx].exposure_time,
            'target_sky_brightness': 0.,
            'target_limiting_mag': 0.,
            'target_metric_value':  0.,
            'target_total_requests_tonight': 1,  
            'request_id': self.queue.index[idx]}

        return next_obs

    def _remove_requests(self, request_id):
        """Remove a request from the queue"""

        if self.queue.loc[request_id,'n_repeats'] > 1:
            self.queue.loc[request_id,'n_repeats'] -= 1
        else:    
            self.queue = self.queue.drop(request_id)

    def _return_queue(self):

        # by construction the list queue is already in order
        queue = self.queue.copy()
        queue['ordered'] = True
        return queue


class RequestPool(object):

    def __init__(self):
        # initialize empty dataframe to add to
        # TODO: currently treating the index as the request_id; should it be
        # unique across sessions?
        self.pool = pd.DataFrame()
        pass

    def add_request_sets(self, program_id, subprogram_name, 
                field_ids, filter_ids, total_requests_tonight, priority=1):
        """program_ids must be scalar"""

        assert (scalar_len(program_id) == 1) 
        assert (scalar_len(subprogram_name) == 1) 

        n_fields = scalar_len(field_ids)
        if n_fields == 1:
            field_ids = [field_ids]

        # build df as a list of dicts
        request_sets = []
        for i, field_id in enumerate(field_ids):
            request_sets.append({
                'program_id': program_id,
                'subprogram_name': subprogram_name,
                'field_id': field_id,
                'filter_ids': filter_ids.copy(),
                'total_requests_tonight': total_requests_tonight,
                'priority': priority})

        self.pool = self.pool.append(pd.DataFrame(request_sets), 
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
        filters.remove(filter_id)
        if len(filters) == 0:
            self.remove_request_sets(request_set_id)
        else:
            self.pool.set_value(request_set_id, 'filter_ids', filters)

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
        ref_obs = df.apply(get_ref_obs_time, args=(current_state,), axis=1)
        dt = current_state['current_time'].mjd - ref_obs
        stats_str += "\t\tMax time to ref_obs for first obs: {}\n".format(
            dt[wfirst].max())
        stats_str += "\t\tMax time to ref_obs for all obs: {}\n".format(
            dt[walt].max())

    return stats_str
