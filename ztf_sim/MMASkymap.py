"""MMA Skymaps."""

import logging
import types

import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time

from .configuration import QueueConfiguration
from .constants import BASE_DIR
from .utils import approx_hours_of_darkness
from .Fields import Fields
from .QueueManager import GreedyQueueManager, RequestPool

class MMASkymap(object):
    """
    A class to handle operations related to Multi-Messenger Astrophysics (MMA) skymaps.

    Attributes:
    -----------
    trigger_name : str
        The name of the trigger event.
    trigger_time : float
        The time of the trigger event.
    skymap_fields : pandas.DataFrame
        A DataFrame containing the skymap fields with 'field_id' and 'probability' columns.
    fields : Fields, optional
        An instance of the Fields class. If not provided, a new instance is created.
    
    Methods:
    --------
    make_queue(validity_window, observing_fraction=0.5):
        Creates an observation queue based on the skymap fields and the specified validity window and observing fraction.
    return_skymap():
        Returns the skymap fields DataFrame.
    persist_skymap():
        Placeholder method to persist the skymap.
    archive_persisted_skymap():
        Placeholder method to archive the persisted skymap.
    """

    def __init__(self, trigger_name, trigger_time, skymap_fields, fields=None):

        self.logger = logging.getLogger(__name__)

        self.trigger_name = trigger_name
        self.trigger_time = trigger_time
        self.skymap_fields = pd.DataFrame(skymap_fields)
        assert('field_id' in self.skymap_fields)
        assert('probability' in self.skymap_fields)

        if fields is None:
            self.fields = Fields()
        else:
            self.fields = fields

    def make_queue(self, validity_window, observing_fraction=0.5):

        assert (0 <= observing_fraction <= 1)

        # use a generic configuration and override
        queue_config = QueueConfiguration(BASE_DIR+'../sims/missed_obs.json')
        queue_name = self.trigger_name+'_greedy'
        queue_config.config['queue_name'] = queue_name
        queue_config.config['queue_description'] = queue_name
        queue_config.config['queue_manager'] = 'greedy'
        queue_config.config['observing_programs'] = []
        queue_config.config['validity_window_mjd'] = validity_window

        Time_validity_start = Time(validity_window[0], format='mjd')

        # visibility check
        self.fields.compute_observability(Time_validity_start)
        observable_field_ids = self.fields.select_field_ids(dec_range=[-32,90.],
                           grid_id=0,
                           # use a minimal observable hours cut
                           observable_hours_range=[0.5, 24.])

        # only select fields that are observable tonight and in the primary grid
        w = self.skymap_fields['field_id'].apply(lambda x: x in observable_field_ids)
        skymap_fields = self.skymap_fields.loc[w,:]

        # sort by probability 
        skymap_field_ids = skymap_fields.sort_values(by='probability', ascending=False)['field_id'].values.tolist()
        
        # limit to # of fields allowed during the night
        # for now we're not going to try to propagate in the exact allocation 
        # of observable time; instead we'll just apply a fraction
        # Let's not assume that the validity range provided is only dark time
        dark_time = approx_hours_of_darkness(Time_validity_start,
                                             twilight=coord.Angle(18*u.degree))

        n_fields = int((dark_time * observing_fraction 
                        / (40*u.second) / 2).to(u.dimensionless_unscaled))

        skymap_field_ids = skymap_field_ids[:n_fields]

        w = skymap_fields['field_id'].apply(lambda x: x in skymap_field_ids)
        skymap_fields = skymap_fields.loc[w,:]


        rp = RequestPool()
        for idx, row in skymap_fields.iterrows():
            rp.add_request_sets(1,
                                'MSIP_EMGW',
                                'Kulkarni',
                                int(row['field_id']),
                                [1,2],
                                30*u.minute,
                                30*u.second,
                                2,
                                probability=row['probability'])

        queue = GreedyQueueManager(queue_name, queue_config, rp = rp)

        self.logger.info(f"""Making queue for {self.trigger_name} with """
                         f"""{[(int(row['field_id']), row['probability']) for idx, row in skymap_fields.iterrows()]}""")

        return queue


    def return_skymap(self):
        return self.skymap_fields

    def persist_skymap(self):
        pass
    
    def archive_persisted_skymap(self):
        pass

