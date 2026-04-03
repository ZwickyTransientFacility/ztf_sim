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
    """Multi-messenger skymap queue builder.

    Wraps a gravitational-wave (or other multi-messenger) probability sky map
    and generates a ``GreedyQueueManager`` follow-up queue sorted by field
    probability.

    Attributes
    ----------
    trigger_name : str
        Unique identifier for the alert event.
    trigger_time : astropy.time.Time
        UTC time of the alert.
    skymap_fields : pandas.DataFrame
        Per-field probabilities with columns ``'field_id'`` and
        ``'probability'``.
    fields : Fields
        ZTF field grid used for observability checks.
    """

    def __init__(self, trigger_name, trigger_time, skymap_fields, fields=None):
        """Initialise an MMA skymap.

        Parameters
        ----------
        trigger_name : str
            Unique identifier for the alert.
        trigger_time : astropy.time.Time
            UTC time of the gravitational-wave (or other) alert.
        skymap_fields : dict or pandas.DataFrame
            Per-field probabilities. Must contain columns ``'field_id'``
            and ``'probability'``.
        fields : Fields or None, optional
            Pre-loaded field grid. If ``None``, a new `Fields` instance is
            created. Default is ``None``.

        Raises
        ------
        AssertionError
            If *skymap_fields* does not contain the required columns.
        """

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
        """Build a greedy follow-up queue for this skymap.

        Selects the highest-probability observable fields up to a number
        proportional to *observing_fraction* of the available dark time,
        sorted by descending probability. Restricts to primary grid (grid_id
        = 0) fields with at least 0.5 hours of observability.

        Parameters
        ----------
        validity_window : list of float
            ``[start_mjd, stop_mjd]`` defining the queue's active window.
        observing_fraction : float, optional
            Fraction of dark time (18-degree twilight) to allocate to this
            event. Default is 0.5.

        Returns
        -------
        GreedyQueueManager
            Follow-up queue containing the selected fields with their skymap
            probabilities as request weights.

        Raises
        ------
        AssertionError
            If *observing_fraction* is not in [0, 1].
        """

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
        """Return the skymap field probability table.

        Returns
        -------
        pandas.DataFrame
            The ``skymap_fields`` DataFrame with columns ``'field_id'`` and
            ``'probability'``.
        """
        return self.skymap_fields

    def persist_skymap(self):
        """Persist the skymap to disk.

        Not yet implemented; reserved for future use.
        """
        pass

    def archive_persisted_skymap(self):
        """Move a persisted skymap to a dated archive directory.

        Not yet implemented; reserved for future use.
        """
        pass

