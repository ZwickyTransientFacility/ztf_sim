
import numpy as np
import pandas as pd
from fields import Fields

class QueueManager:

    def __init__(self, rp = None, fields = None):
        if rp is None:
            # initialize an empty RequestPool
            self.rp = RequestPool()
        else:
            self.rp = rp

        if fields is None:
            self.fields = Fields()
        else:
            self.fields = fields

    def next_obs(self, current_state):
        """Given current state, return the parameters for the next request"""
        # don't store the telescope state locally!
        
        # define functions that actually do the work in subclasses
        return self._next_obs(current_state)

    def update_queue(self, **kwargs):
        """Recalculate queue"""

        # define functions that actually do the work in subclasses
        return self._update_queue()


class GreedyQueueManager(QueueManager):

    def __init__(self, **kwargs):
        super(GreedyQueueManager, self).__init__(**kwargs)


    def _next_obs(self, current_state):
        """Select the highest value request."""

        # since this is a greedy queue, we update the queue after each obs
        self._update_queue(current_state)

        # request_id of the highest value request
        max_idx = self.queue.value.argmax()
        
        return {'target_field_id': self.queue.ix[max_idx].field_id,
                'target_ra': self.queue.ix[max_idx].ra,
                'target_dec': self.queue.ix[max_idx].dec,
                'target_filter_id': self.queue.ix[max_idx].filter_id,
                'target_program_id': self.queue.ix[max_idx].program_id,
                'request_id': max_idx}

    def _update_queue(self, current_state):
        """Calculate greedy weighting of requests in the Pool using current 
        telescope state only"""

        # join with fields so we have the information we need
        df = self.rp.pool.join(self.fields.fields,on='field_id')

        # initialize cadence windows
        df['in_cadence_window'] = False

        # use cadence functions to compute requests with active cadence windows
        for idx, row in df.iterrows():
            df['in_cadence_window'].ix[idx] = \
                    eval('{}(row, current_state)'.format(row['cadence_func']))

        cadence_cuts = df['in_cadence_window'] 
        df = df[cadence_cuts]

        # compute readout/slew overhead times, plus current altitude
        # TODO: need to add overhead for filter changes
        df_overhead, df_alt = self.fields.overhead_time(current_state, 
                cuts=cadence_cuts)

        df = df.join(df_overhead, left_on='field_id')

        # airmass cut (or add airmass weighting to value below)
        airmass_cuts = zenith_angle_to_airmass(90.*u.deg - df_alt) <= MAX_AIRMASS 
        df = df[airmass_cuts]

        # TODO: penalize volume for both extinction (airmass) and fwhm penalty
        # due to atmospheric refraction (use a lookup table, probably)
        # as well as moon phase and distance
        
        # value = 1./(obs_time + slew time)
        # TODO: for now, make this unitless
        df['value'] = 1.*u.sec/(EXPOSURE_TIME + df['overhead_time'])
        
        self.queue = df


class RequestPool:

    def __init__(self):
        # initialize empty dataframe to add to
        # TODO: currently treating the index as the request_id; should it be
        # unique across sessions?
        self.pool = pd.DataFrame()
        pass

    def add_requests(self, program_id, field_ids, filter_id, 
            cadence_func, cadence_pars, priority=1):
        """all scalars except field_ids"""
        #TODO: Compound Requests

        def scalar_len(x):
            return len(np.atleast_1d(x))

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
                'priority': priority})

        self.pool = self.pool.append(pd.DataFrame(requests),ignore_index=True)

    def n_requests(self):
        return len(self.pool)

    def remove_requests(self, request_ids):
        """Remove completed or otherwise unwanted requests by request_id
        
        request_ids : scalar or list
            requests to drop (index of self.pool)"""
        self.pool = self.pool.drop(request_ids)
