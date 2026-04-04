"""Classes implementing Observing Programs."""

import logging
import numpy as np
import astropy.units as u
from astropy.time import Time
from .constants import EXPOSURE_TIME, READOUT_TIME, TIME_BLOCK_SIZE
from .utils import approx_hours_of_darkness
from .field_selection_functions import *


class ObservingProgram(object):
    """Encodes the field selection, cadence, and time allocation for one science program.

    Called nightly by `QueueManager.assign_nightly_requests` to generate the
    set of observation requests for tonight. Field selection may be static (an
    explicit list of field IDs or a position-based cut) or dynamic (a named
    function from ``field_selection_functions.py``).

    Attributes
    ----------
    program_id : int
        Program identifier (0 = engineering, 1 = MSIP, 2 = collaboration,
        3 = Caltech).
    subprogram_name : str
        Human-readable subprogram label, unique within the program.
    filter_ids : list of int
        Filter sequence for observations (1 = g, 2 = r, 3 = i).
    internight_gap : astropy.units.Quantity
        Minimum time between successive nightly observations of the same field.
    intranight_gap : astropy.units.Quantity
        Minimum time between repeat visits to the same field within one night.
    n_visits_per_night : int
        Number of visits per field per night.
    """

    def __init__(self, program_id, subprogram_name, program_pi,
                 program_observing_time_fraction, subprogram_fraction,
                 field_ids, filter_ids, internight_gap,
                 intranight_gap, n_visits_per_night,
                 exposure_time = EXPOSURE_TIME,
                 nobs_range=None,
                 filter_choice='rotate',
                 active_months='all',
                 field_selection_function=None):
        """Initialise an observing program.

        Exactly one of *field_ids* and *field_selection_function* must be
        provided (not both, not neither).

        Parameters
        ----------
        program_id : int
            Program identifier (0–3, see ``PROGRAM_NAME_TO_ID``).
        subprogram_name : str
            Human-readable subprogram label.
        program_pi : str
            Principal investigator name.
        program_observing_time_fraction : float
            Fraction of total telescope time allocated to this program (0–1).
        subprogram_fraction : float
            Fraction of the program's time allocated to this subprogram (0–1).
        field_ids : list of int or None
            Static list of ZTF field IDs. Mutually exclusive with
            *field_selection_function*.
        filter_ids : list of int
            Filters to observe, e.g. ``[1, 2]``. Interpretation depends on
            *filter_choice*.
        internight_gap : astropy.units.Quantity
            Minimum elapsed time between observations of the same field on
            successive nights.
        intranight_gap : astropy.units.Quantity
            Minimum elapsed time between repeat visits to the same field
            within one night.
        n_visits_per_night : int
            Number of visits per field per night.
        exposure_time : astropy.units.Quantity, optional
            Per-visit exposure time. Default is ``EXPOSURE_TIME`` (30 s).
        nobs_range : dict or None, optional
            If provided, restricts eligible fields to those whose observation
            count falls within ``{'min_obs': int, 'max_obs': int}``.
            Additional optional keys: ``'program_ids'``, ``'subprogram_names'``,
            ``'filter_ids'``, ``'mjd_range'``.
        filter_choice : str, optional
            ``'rotate'``: use one filter per night, cycling through
            *filter_ids* by ``floor(mjd) % n_filters``.
            ``'sequence'``: use the complete *filter_ids* sequence every
            night (length must equal *n_visits_per_night*).
            Default is ``'rotate'``.
        active_months : str or list of int, optional
            ``'all'`` or a list of month numbers (1–12) during which this
            program is active. Default is ``'all'``.
        field_selection_function : str or None, optional
            Name of a function in ``field_selection_functions.py`` that
            dynamically selects fields each night. The special value
            ``'EP-bypass'`` signals that this program is handled by
            ``make_nightly_timed_blocks`` and should return an empty list.
            Mutually exclusive with *field_ids*.

        Raises
        ------
        AssertionError
            If both or neither of *field_ids* and *field_selection_function*
            are provided.
        """

        assert ((field_ids is None) or (field_selection_function is None))
        assert not((field_ids is None) and (field_selection_function is None))

        self.logger = logging.getLogger(__name__)

        self.program_id = program_id
        self.subprogram_name = subprogram_name
        self.program_pi = program_pi
        self.program_observing_time_fraction = program_observing_time_fraction
        self.subprogram_fraction = subprogram_fraction
        self.field_ids = field_ids
        self.filter_ids = filter_ids

        self.internight_gap = internight_gap
        self.intranight_gap = intranight_gap
        self.n_visits_per_night = n_visits_per_night
        self.exposure_time = exposure_time # a Quantity

        self.nobs_range = nobs_range 
        self.filter_choice = filter_choice
        if active_months != 'all':
            # allow scalar or list input
            self.active_months = np.atleast_1d(active_months)
        else:
            self.active_months = 'all'

        self.field_selection_function = field_selection_function

    def assign_nightly_requests(self, time, fields, obs_log,
            other_program_fields,
            block_programs=False, skymaps = None, **kwargs):
        """Generate tonight's observation request sets for this program.

        Applies a five-step pipeline:

        1. **Month guard** — returns ``[]`` if the current month is not in
           ``active_months``.
        2. **Observable fields** — finds fields visible for at least
           ``n_visits_per_night`` consecutive 30-min blocks.
        3. **Cadence filter** — if using a static field list, removes fields
           observed within ``internight_gap − 0.6 day``; if using a dynamic
           selection function, cadence is handled internally by that function.
        4. **nobs_range filter** — if ``nobs_range`` is configured, retains
           only fields with ``min_obs ≤ n_obs ≤ max_obs``.
        5. **Filter sequence** — constructs the per-visit filter list for
           tonight using ``filter_choice``.

        Parameters
        ----------
        time : astropy.time.Time
            Start of the current night (used for month/cadence checks).
        fields : Fields
            ZTF field grid object (must have ``compute_blocks`` and
            ``compute_observability`` called for tonight).
        obs_log : ObsLogger
            Observation history for cadence queries.
        other_program_fields : dict
            Mapping ``(program_id, subprogram_name) -> dict`` with keys
            ``'field_ids'``, ``'field_selection_function'``, and
            ``'requests_allowed'``. Used by dynamic selection functions to
            avoid field overlap between programs.
        block_programs : bool, optional
            Reserved for future use. Default is ``False``.
        skymaps : dict or None, optional
            Mapping of skymap name to skymap object, passed to dynamic
            selection functions.
        **kwargs
            Passed to the dynamic field selection function if one is used.

        Returns
        -------
        list of dict
            A one-element list containing a request-set dict with keys:
            ``'program_id'``, ``'subprogram_name'``, ``'program_pi'``,
            ``'field_ids'`` (numpy.ndarray), ``'filter_ids'`` (list),
            ``'exposure_time'`` (Quantity), ``'intranight_gap'`` (Quantity),
            ``'total_requests_tonight'`` (int). Returns ``[]`` if no fields
            are eligible tonight.
        """

        # filters are given in filter_ids:
        # either a set of filters, or a fixed sequence
        # filter_choice = 'rotate':
        #   use one filter set per night, keyed by mjd % n_filters
        # filter_choice = 'sequence':
        #   use hard-coded sequence given in filter_ids

        # determine if this program is active this month and return
        # an empty set if not
        if self.active_months != 'all':
            if time.to_datetime().month not in self.active_months:
                return []
        
        
        # compute nightly altaz blocks and observability windows
        fields.compute_blocks(time)
        fields.compute_observability(time)

        n_filters = len(set(self.filter_ids))
        if self.filter_choice == 'rotate':
            night_index_filters = np.floor(time.mjd % n_filters).astype(int)
            filter_ids_tonight = self.filter_ids[night_index_filters]
            filter_ids_last_night = self.filter_ids[night_index_filters - 1]
            # make it a list
            filter_ids_tonight = [filter_ids_tonight]
        else:
            filter_ids_tonight = list(set(self.filter_ids))

        # Choose which fields will be observed

        # minimum time to observe N visits
        obs_field_ids = fields.select_field_ids(observable_hours_range=
            [(self.n_visits_per_night*TIME_BLOCK_SIZE).to(u.hour).value, 24.])

        # if needed, compute the OP fields on a nightly basis
        if self.field_selection_function is not None:
            if self.field_selection_function == 'EP-bypass':
                return []
            try:
                selection_function = globals()[self.field_selection_function]
                field_ids = selection_function(time, obs_log, other_program_fields, fields, skymaps)
                self.logger.info(f'Program ID {self.program_id}, subprogram {self.subprogram_name}: selected {len(field_ids)} fields')
                self.logger.debug(f'    {field_ids}')
            except Exception as e:
                #raise(e) # needed to debug filter_selection_functions
                self.logger.exception(e)
                self.logger.warning(f'Error in generating nightly field list for Program ID {self.program_id}, subprogram {self.subprogram_name}, returning zero fields!')  
                return []
        else:
            field_ids = self.field_ids

        # now form the intersection of observable fields and the OP fields
        pool_ids = obs_field_ids.intersection(field_ids)
        self.logger.debug(f'Program ID {self.program_id}, subprogram {self.subprogram_name}: {len(pool_ids)} fields observable')

        # get the times they were last observed:
        # (note that fields *never* observed will not be included)
        # since this is function is for determining requests
        # at the start of the night, exclude observations taken tonight
        # this lets us restart the scheduler without breaking things
        last_observed_times = obs_log.select_last_observed_time_by_field(
                field_ids = pool_ids,
                filter_ids = filter_ids_tonight,
                program_ids = [self.program_id],
                subprogram_names = [self.subprogram_name],
                # arbitrary early date; start of night tonight
                mjd_range = [Time('2001-01-01').mjd,np.floor(time.mjd)])

        # we want an object observed at the end of the night N days ago
        # to be observed at the start of the night now.
        # Max night length is 12.2 hours
        cutoff_time = (time - (self.internight_gap - 0.6 * u.day)).mjd

        # find fields last observed more recently than that
        wrecent = (last_observed_times['expMJD'] >= cutoff_time)
        recent_field_ids = last_observed_times.loc[wrecent].index.tolist()

        if self.field_selection_function is None:
            # reduce the list to only those not recently observed:
            pool_ids_old = [idi for idi in pool_ids if idi not in recent_field_ids]
            request_fields = fields.fields.loc[pool_ids_old]
        else:
            # field_selection_function needs to apply the cadence cut
            request_fields = fields.fields.loc[pool_ids]

        # if we have an nobs_range argument (eg for reference building), use it
        if self.nobs_range is not None:
            if 'program_ids' not in self.nobs_range:
                program_ids = None
            else:
                program_ids = self.nobs_range['program_ids']

            if 'subprogram_names' not in self.nobs_range:
                subprogram_names = None
            else:
                subprogram_names = self.nobs_range['subprogram_names']

            if 'filter_id' in self.nobs_range:
                self.nobs_range['filter_ids'] = [self.nobs_range['filter_id']]
            if 'filter_ids' not in self.nobs_range:
                filter_ids = None
            else:
                filter_ids = self.nobs_range['filter_ids'] 
                
            if 'mjd_range' not in self.nobs_range:
                mjd_range = None
            else:
                mjd_range = self.nobs_range['mjd_range']

            assert 'min_obs' in self.nobs_range
            assert 'max_obs' in self.nobs_range
                
            nobs = obs_log.select_n_obs_by_field(filter_ids = filter_ids,
                    program_ids = program_ids, 
                    subprogram_names = subprogram_names,
                    mjd_range = mjd_range)
            
            # function above only returns fields that have been observed at
            # least once.  use the intersection if min_obs > 0:
            w = ((nobs >= self.nobs_range['min_obs']) & 
                    (nobs <= self.nobs_range['max_obs']))
            if self.nobs_range['min_obs'] > 0:
                nobs_inrange = nobs.loc[w]
                request_fields = request_fields.join(nobs_inrange,how='inner')
            else:
                # drop rows out of range (which here means only those with 
                # nobs > max_obs
                nobs_outofrange = nobs.loc[~w]
                # find fields that are in request_fields but out of range
                nobs_outofrange = request_fields.join(nobs_outofrange,how='inner')
                # now drop them
                request_fields = request_fields.drop(nobs_outofrange.index)
            

        # construct request sets: list of inputs to RequestPool.add_requests
        # scalar everything except field_ids

        if self.filter_choice == 'rotate':
            filter_sequence = [filter_ids_tonight[0] for i in
                               range(self.n_visits_per_night)]
        elif self.filter_choice == 'sequence':
            assert(len(self.filter_ids) == self.n_visits_per_night)
            filter_sequence = self.filter_ids.copy()

        self.logger.debug(f'Program ID {self.program_id}, subprogram {self.subprogram_name}: {len(request_fields.index.values)} fields requested')

        request_set = []
        request_set.append(
            {'program_id': self.program_id,
             'subprogram_name': self.subprogram_name,
             'program_pi': self.program_pi,
             'field_ids': request_fields.index.values,
             'filter_ids': filter_sequence,
             'exposure_time': self.exposure_time,
             'intranight_gap': self.intranight_gap,
             'total_requests_tonight': self.n_visits_per_night})

        return request_set

    def time_per_exposure(self):
        """Return the total time consumed per exposure including readout.

        Returns
        -------
        astropy.units.Quantity
            ``exposure_time + READOUT_TIME``.
        """
        return self.exposure_time + READOUT_TIME
