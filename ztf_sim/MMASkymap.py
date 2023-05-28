"""MMA Skymaps."""

import logging
import types

import pandas as pd
import astropy.units as u

from .configuration import QueueConfiguration
from .QueueManager import GreedyQueueManager, RequestPool
from .constants import BASE_DIR

class MMASkymap(object):

    def __init__(self, trigger_name, trigger_time, skymap_fields, fields=None):

        self.logger = logging.getLogger(__name__)

        self.trigger_name = trigger_name
        self.trigger_time = trigger_time
        self.skymap_fields = pd.DataFrame(skymap_fields)
        assert('field_id' in self.skymap_fields)
        assert('probability' in self.skymap_fields)

#        if fields is None:
#            self.fields = Fields()
#        else:
#            self.fields = fields

    def make_queue(self, validity_window):

        # use a generic configuration and override
        queue_config = QueueConfiguration(BASE_DIR+'../sims/missed_obs.json')
        queue_name = self.trigger_name+'_greedy'
        queue_config.config['queue_name'] = queue_name
        queue_config.config['queue_description'] = queue_name
        queue_config.config['queue_manager'] = 'greedy'
        queue_config.config['observing_programs'] = []
        queue_config.config['validity_window_mjd'] = validity_window

        rp = RequestPool()
        for idx, row in self.skymap_fields.iterrows():
            rp.add_request_sets(1,
                                'MSIP_EMGW',
                                'Kulkarni',
                                row['field_id'],
                                [1,2],
                                30*u.minute,
                                30*u.second,
                                2,
                                probability=row['probability'])

        queue = GreedyQueueManager(queue_name, queue_config, rp = rp)

        return queue


    def return_skymap(self):
        return self.skymap_fields

    def persist_skymap(self):
        pass
    
    def archive_persisted_skymap(self):
        pass

