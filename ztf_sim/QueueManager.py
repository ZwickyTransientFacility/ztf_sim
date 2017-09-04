from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import object
import numpy as np
import pandas as pd
from fields import Fields
from sky_brightness import SkyBrightness, FakeSkyBrightness
from magnitudes import limiting_mag
import astropy.coordinates as coord
from cadence import *
from optimize import request_set_optimize, slot_optimize, tsp_optimize
from constants import *
from utils import *
import pdb

class QueueEmptyError(Exception):
    """Error class for when the nightly queue has no more fields"""
    pass


class QueueManager(object):

    def __init__(self, observing_programs=[], rp=None, fields=None,
                 block_programs=True):

        # list of ObservingPrograms
        self.observing_programs = observing_programs

        # block on which the queue parameters were calculated
        self.queue_slot = None

        # number allowed requests by program tonight (dict)
        self.requests_allowed = {}

        # count of executed requests by program tonight (dict)
        self.requests_completed = {id:0 for id in PROGRAM_IDS}

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

    def assign_nightly_requests(self, current_state):
        # clear previous request pool
        self.rp.clear_all_request_sets()
        # reset the first observation of the night counters
        self.fields.clear_first_obs()
        # clear count of executed requests 
        self.requests_completed = {id:0 for id in PROGRAM_IDS}
        # set number of allowed requests by program.
        self.determine_allowed_requests(current_state['current_time'])

        for program in self.observing_programs:

            request_sets = program.assign_nightly_requests(
                current_state['current_time'], self.fields,
                block_programs=self.block_programs)
            for rs in request_sets:
                self.rp.add_request_sets(rs['program_id'], rs['field_ids'],
                                     rs['filter_ids'], 
                                     rs['total_requests_tonight'])

        assert(len(self.rp.pool) > 0)

        # any specific tasks needed)
        self._assign_nightly_requests(current_state)

    def determine_allowed_requests(self, time):
        """Use count of past observations and expected observing time fractions
        to determine number of allowed requests tonight."""

        self.requests_allowed = {id:0 for id in PROGRAM_IDS}
        
        obs_count_by_program = self.fields.count_total_obs_by_program()
        total_obs = np.sum(list(obs_count_by_program.values()))
        for program in self.observing_programs:
            n_requests = program.number_of_allowed_requests(time)
            delta = np.round(
                (obs_count_by_program[program.program_id] -
                 program.program_observing_time_fraction * total_obs)
                 * program.subprogram_fraction)

            # TODO: tweak as needed
            # how quickly do we want to take to reach equalization?
            CATCHUP_FACTOR = 0.20
            n_requests -= np.round(delta * CATCHUP_FACTOR).astype(np.int)
            self.requests_allowed[program.program_id] += n_requests

        for id, n_requests in list(self.requests_allowed.items()):
            if n_requests < 0:
                self.requests_allowed[id] = 0

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
        """Remove a request from both the queue and the request pool"""

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
        df.loc[:, 'moonillf'] = astroplan.moon.moon_illumination(
            # Don't use P48_loc to avoid astropy bug:
            # https://github.com/astropy/astroplan/pull/213
            time)
        # time, P48_loc)
        df.loc[:, 'moon_dist'] = sc.separation(moon).to(u.deg).value
        df.loc[:, 'moonalt'] = moon_altaz.alt.to(u.deg).value
        df.loc[:, 'sunalt'] = sun_altaz.alt.to(u.deg).value

        # compute sky brightness
        # only have values for reasonable altitudes (set by R20_absorbed...)
        wup = df['altitude'] > 10
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
        # assign a very bright limiting mag to the fields that are down 
        # so the metric goes to zero
        df.loc[~wup, 'limiting_mag'] = -99

        return df['limiting_mag']



class GurobiQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GurobiQueueManager, self).__init__(**kwargs)
        self.block_obs_number = 0

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

        # TODO: make the queue have the right datatypes
        next_obs = {'target_field_id': int(self.queue.ix[idx].field_id),
            'target_ra': self.queue.ix[idx].ra,
            'target_dec': self.queue.ix[idx].dec,
            'target_filter_id': int(self.filter_by_slot[self.queue_slot]),
            'target_program_id': int(self.queue.ix[idx].program_id),
            'target_exposure_time': EXPOSURE_TIME,
            'target_sky_brightness': 0.,
            'target_limiting_mag': 0.,
            'target_metric_value':  0.,
            'target_total_requests_tonight':
            int(self.queue.ix[idx].total_requests_tonight),
            'request_id': idx}

#            'target_sky_brightness': self.queue.ix[idx].sky_brightness,
#            'target_limiting_mag': self.queue.ix[idx].limiting_mag,
#            'target_metric_value':  self.queue.ix[idx].value,
#            'target_request_number_tonight':
#            'target_total_requests_tonight':
#            self.queue.ix[idx].total_requests_tonight,
#            'request_id': idx}

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
                lim_mags[(bi, fid)] = \
                    self.compute_limiting_mag(df, ti, filter_id = fid)

        # this results in a MultiIndex on the *columns*: level 0 is block,
        # level 1 is filter_id.  df_metric.unstack() flattens it
        self.block_lim_mags = pd.DataFrame(lim_mags)
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

        self.queue_slot = block_index(current_state['current_time'])

        # retrieve requests to be observed in this block
        req_list = self.queued_requests_by_slot[self.queue_slot].values[0]

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
            df[self.queue_slot], axis='dome')
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

    def _remove_requests(self, request_id):
        """Remove a request from both the queue and the pool"""

        # should be the topmost item
        assert (self.queue_order[0] == request_id)
        self.queue_order = self.queue_order[1:]
        self.queue = self.queue.drop(request_id)
        # TODO: decrement it in the pool


class GreedyQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GreedyQueueManager, self).__init__(**kwargs)

    def _assign_nightly_requests(self, current_state):
        pass

    def _next_obs(self, current_state):
        """Select the highest value request."""

        # since this is a greedy queue, we update the queue after each obs
        # for speed, only do the whole recalculation if we're in a new slot
        if ((block_index(current_state['current_time'])[0] != self.queue_slot)
                or (len(self.queue) == 0)):
            self._update_queue(current_state)
        else:
            # otherwise just recalculate the overhead times
            _ = self._update_overhead(current_state)

        # request_id of the highest value request
        max_idx = self.queue.value.argmax()

        # enforce program balance
        progid = self.queue.ix[max_idx].program_id
        if ((self.requests_completed[progid] + 1) >= 
                self.requests_allowed[progid]):
            # this program has used up all its obs for tonight.
            # remove all requests for this program from the pool 
            w = self.rp.pool.program_id == progid
            self.rp.remove_requests(self.rp.pool[w].index.tolist())
            # reset the queue
            self._update_queue(current_state)
            # and request a new next_obs
            next_obs = self._next_obs(current_state)
        else:
            next_obs = {'target_field_id': self.queue.ix[max_idx].field_id,
                'target_ra': self.queue.ix[max_idx].ra,
                'target_dec': self.queue.ix[max_idx].dec,
                'target_filter_id': self.queue.ix[max_idx].filter_id,
                'target_program_id': self.queue.ix[max_idx].program_id,
                'target_exposure_time': EXPOSURE_TIME,
                'target_sky_brightness': self.queue.ix[max_idx].sky_brightness,
                'target_limiting_mag': self.queue.ix[max_idx].limiting_mag,
                'target_metric_value':  self.queue.ix[max_idx].value,
                'target_total_requests_tonight':
                self.queue.ix[max_idx].total_requests_tonight,
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
        self.queue_slot = block_index(current_state['current_time'])

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
        in_window = {}
        for idx, row in df.iterrows():
            # this is way, way slower
            # df['in_cadence_window'].ix[idx] = \
            in_window[idx] = eval(
                '{}(row, current_state)'.format(row['cadence_func']))

        cadence_cuts = pd.Series(in_window)
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

        # compute inputs for sky brightness
        sc = coord.SkyCoord(df['ra'], df['dec'], frame='icrs', unit='deg')
        sun = coord.get_sun(current_state['current_time'])
        sun_altaz = skycoord_to_altaz(sun, current_state['current_time'])
        moon = coord.get_moon(current_state['current_time'],
                              location=P48_loc)
        moon_altaz = skycoord_to_altaz(moon, current_state['current_time'])
        df.loc[:, 'moonillf'] = astroplan.moon.moon_illumination(
            # Don't use P48_loc to avoid astropy bug:
            # https://github.com/astropy/astroplan/pull/213
            current_state['current_time'])
        # current_state['current_time'], P48_loc)
        df.loc[:, 'moon_dist'] = sc.separation(moon).to(u.deg).value
        df.loc[:, 'moonalt'] = moon_altaz.alt.to(u.deg).value
        df.loc[:, 'sunalt'] = sun_altaz.alt.to(u.deg).value

        # compute sky brightness
        df.loc[:, 'sky_brightness'] = self.Sky.predict(df)
        #df = pd.merge(df, df_sky, left_on='field_id', right_index=True)

        # compute seeing at each pointing
        df.loc[:, 'seeing'] = seeing_at_pointing(current_state['current_zenith_seeing'],
                                                 df['altitude'])
        #df_seeing.name = 'seeing'
        #df = pd.merge(df, df_seeing, left_on='field_id', right_index=True)

        df.loc[:, 'limiting_mag'] = limiting_mag(EXPOSURE_TIME, df['seeing'],
                                                 df['sky_brightness'],
                                                 filter_id=df['filter_id'],
                                                 altitude=df['altitude'], SNR=5.)
        #df_limmag.name = 'limiting_mag'
        #df = pd.merge(df, df_limmag, left_on='field_id', right_index=True)

        df.loc[:, 'value'] = self._metric(df)

        self.queue = df

    def _remove_requests(self, request_id):
        """Remove a request from both the queue and the request pool"""

        self.queue = self.queue.drop(request_id)
        raise NotImplementedError("Need to adjust pool behavior to remove single obs from request_sets")
        #self.rp.remove_requests(request_id)



class ListQueueManager(QueueManager):
    """Simple Queue that returns observations in order."""

    def __init__(self, **kwargs):
        super(ListQueueManager, self).__init__(**kwargs)

    def _assign_nightly_requests(self, current_state):
        raise NotImplementedError("ListQueueManager should be loaded by load_queue()")

    def load_queue(self, queue_dict_list):
        """Initialize an ordered queue.

        queue_dict_list is a list of dicts, one per observation"""
        
        self.queue = pd.DataFrame(queue_dict_list)

        # check that major columns are included
        required_columns = ['field_id','ra','dec', 'program_id', 'filter_id']

    def _next_obs(self, current_state):
        """Return the next observation in the time ordered queue unless it has expired."""

        
        if len(self.queue) == 0:
            raise QueueEmptyError("No more observations in queue!")
        
        # take the next observation in line
        idx = 0

        # TODO: check if it has expired
       
        # if no exposure time is specified, use the default
        if 'exposure_time' not in self.queue.columns:
            texp = EXPOSURE_TIME
        else:
            texp =  self.queue.iloc[idx].exposure_time

        next_obs = {'target_field_id': int(self.queue.iloc[idx].field_id),
            'target_ra': self.queue.iloc[idx].ra,
            'target_dec': self.queue.iloc[idx].dec,
            'target_filter_id': self.queue.iloc[idx].filter_id,
            'target_program_id': int(self.queue.iloc[idx].program_id),
            'target_exposure_time': texp,
            'target_sky_brightness': 0.,
            'target_limiting_mag': 0.,
            'target_metric_value':  0.,
            'target_total_requests_tonight': 1,  
            'request_id': self.queue.index[0]}

        return next_obs

    def _remove_requests(self, request_id):
        """Remove a request from the queue"""

        self.queue = self.queue.drop(request_id)


class RequestPool(object):

    def __init__(self):
        # initialize empty dataframe to add to
        # TODO: currently treating the index as the request_id; should it be
        # unique across sessions?
        self.pool = pd.DataFrame()
        pass

    def add_request_sets(self, program_id, field_ids, filter_ids,
                     total_requests_tonight, 
                     priority=1):
        """program_ids must be scalar"""

        assert (scalar_len(program_id) == 1) 

        n_fields = scalar_len(field_ids)
        if n_fields == 1:
            field_ids = [field_ids]

        # build df as a list of dicts
        request_sets = []
        for i, field_id in enumerate(field_ids):
            request_sets.append({
                'program_id': program_id,
                'field_id': field_id,
                'filter_ids': filter_ids,
                'total_requests_tonight': total_requests_tonight,
                'priority': priority})

        self.pool = self.pool.append(pd.DataFrame(request_sets), 
            ignore_index=True)

    def n_request_sets(self):
        return len(self.pool)

    def remove_request_sets(self, request_ids):
        """Remove completed or otherwise unwanted requests by request_id

        request_ids : scalar or list
            requests to drop (index of self.pool)"""
        self.pool = self.pool.drop(request_ids)

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
