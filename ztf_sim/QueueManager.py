
import numpy as np
import pandas as pd
from fields import Fields
from sky_brightness import SkyBrightness, FakeSkyBrightness
from magnitudes import limiting_mag
import astropy.coordinates as coord
from cadence import *
from constants import *
from utils import *


class QueueEmptyError(Exception):
    """Error class for when the nightly queue has no more fields"""
    pass


class QueueManager(object):

    def __init__(self, observing_programs=[], rp=None, fields=None,
                 block_programs=True):

        # list of ObservingPrograms
        self.observing_programs = observing_programs

        # block on which the queue parameters were calculated
        self.queue_block = None

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
        self.rp.clear_all_requests()
        # reset the first observation of the night counters
        self.fields.clear_first_obs()

        for program in self.observing_programs:

            request_sets = program.assign_nightly_requests(
                current_state['current_time'], self.fields,
                block_programs=self.block_programs)
            for rs in request_sets:
                self.rp.add_requests(rs['program_id'], rs['field_ids'],
                                     rs['filter_id'], rs['cadence_func'],
                                     rs['cadence_pars'],
                                     rs['request_number_tonight'],
                                     rs['total_requests_tonight'],
                                     priority=rs['priority'])

    def next_obs(self, current_state):
        """Given current state, return the parameters for the next request"""
        # don't store the telescope state locally!

        # define functions that actually do the work in subclasses
        return self._next_obs(current_state)

    def update_queue(self, **kwargs):
        """Recalculate queue"""

        # define functions that actually do the work in subclasses
        return self._update_queue()

    def remove_requests(self, request_id):
        """Remove a request from both the queue and the request pool"""

        self.queue.drop(request_id, inplace=True)
        self.rp.remove_requests(request_id)


class GreedyQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GreedyQueueManager, self).__init__(**kwargs)

    def _next_obs(self, current_state):
        """Select the highest value request."""

        # since this is a greedy queue, we update the queue after each obs
        # for speed, only do the whole recalculation if we're in a new block
        if ((block_index(current_state['current_time'])[0] != self.queue_block)
                or (len(self.queue) == 0)):
            self._update_queue(current_state)
        else:
            # otherwise just recalculate the overhead times
            _ = self._update_overhead(current_state)

        # request_id of the highest value request
        max_idx = self.queue.value.argmax()

        return {'target_field_id': self.queue.ix[max_idx].field_id,
                'target_ra': self.queue.ix[max_idx].ra,
                'target_dec': self.queue.ix[max_idx].dec,
                'target_filter_id': self.queue.ix[max_idx].filter_id,
                'target_program_id': self.queue.ix[max_idx].program_id,
                'target_exposure_time': EXPOSURE_TIME,
                'target_sky_brightness': self.queue.ix[max_idx].sky_brightness,
                'target_limiting_mag': self.queue.ix[max_idx].limiting_mag,
                'target_metric_value':  self.queue.ix[max_idx].value,
                'target_request_number_tonight':
                self.queue.ix[max_idx].request_number_tonight,
                'target_total_requests_tonight':
                self.queue.ix[max_idx].total_requests_tonight,
                'request_id': max_idx}

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
        self.queue_block = block_index(current_state['current_time'])

        # join with fields so we have the information we need
        df = self.rp.pool.join(self.fields.fields, on='field_id')

        df = self._update_overhead(current_state, df=df)

        # start with conservative altitude cut;
        # airmass weighting applied naturally below
        df = df[df['altitude'] > 20]

        if len(df) == 0:
            raise QueueEmptyError("No fields in queue")

        # if restricting to one program per block, drop other programs
        if self.block_programs:
            current_block_program = PROGRAM_BLOCK_SEQUENCE[
                self.queue_block % LEN_BLOCK_SEQUENCE]
            df = df[df['program_id'] == current_block_program]

        # use cadence functions to compute requests with active cadence windows
        # this is slow, so do it after our altitude cut
        in_window = {}
        for idx, row in df.iterrows():
            # this is way, way slower
            # df['in_cadence_window'].ix[idx] = \
            in_window[idx] = eval(
                '{}(row, current_state)'.format(row['cadence_func']))

        cadence_cuts = pd.Series(in_window)
        # TODO: handle if cadence cuts returns no fields
        if np.sum(cadence_cuts) == 0:
            raise QueueEmptyError("No fields with observable cadence windows")
        df = df[cadence_cuts]

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


class RequestPool(object):

    def __init__(self):
        # initialize empty dataframe to add to
        # TODO: currently treating the index as the request_id; should it be
        # unique across sessions?
        self.pool = pd.DataFrame()
        pass

    def add_requests(self, program_id, field_ids, filter_id,
                     cadence_func, cadence_pars, request_number_tonight,
                     total_requests_tonight, priority=1):
        """all scalars except field_ids"""
        # TODO: Compound Requests

        assert ((scalar_len(program_id) == 1) and
                (scalar_len(filter_id) == 1) and
                (scalar_len(cadence_func) == 1))

        n_fields = scalar_len(field_ids)
        if n_fields == 1:
            field_ids = [field_ids]

        # build df as a list of dicts
        requests = []
        for i, field_id in enumerate(field_ids):
            requests.append({
                'program_id': program_id,
                'field_id': field_id,
                'filter_id': filter_id,
                'cadence_func': cadence_func,
                'cadence_pars': cadence_pars,
                'request_number_tonight': request_number_tonight,
                'total_requests_tonight': total_requests_tonight,
                'priority': priority})

        self.pool = self.pool.append(pd.DataFrame(requests), ignore_index=True)

    def n_requests(self):
        return len(self.pool)

    def remove_requests(self, request_ids):
        """Remove completed or otherwise unwanted requests by request_id

        request_ids : scalar or list
            requests to drop (index of self.pool)"""
        self.pool = self.pool.drop(request_ids)

    def clear_all_requests(self):
        self.pool = pd.DataFrame()
