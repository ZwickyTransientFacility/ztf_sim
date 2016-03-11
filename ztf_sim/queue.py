
import pandas as pd

class QueueManager:

    def __init__(self, rp = None):
        if rp is None:
            # initialize an empty RequestPool
            self.rp = RequestPool()
        else:
            self.rp = rp

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
        self.pool = pd.DataFrame()
        pass

    def add_requests(self, program_id, field_id, filter_id, 
            cadence_func, cadence_pars, priority=1):
        #TODO: Compound Requests
        pass

    def remove_requests(self, request_ids):
        """Remove completed or otherwise unwanted requests by request_id"""
        pass
