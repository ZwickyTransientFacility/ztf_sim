
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
        
        return {'target_fieldid': fieldid,
                'target_ra': ra,
                'target_dec': dec,
                'target_filterid': filterid,
                'target_program': program,
                'request_id': request_id}

    def _update_queue(self):
        """Calculate greedy weighting of requests in the Pool using current 
        telescope state only"""

        # initialize all to zero so fields below the horizon aren't counted
        self.rp.pool['value'] = 0.

        # select requests with active cadence windows

        # zenith cut (or airmass weighting)
        
        # slew time
        pass


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
