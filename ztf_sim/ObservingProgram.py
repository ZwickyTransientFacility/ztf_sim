from constants import *
from utils import approx_hours_of_darkness


class ObservingProgram(object):

    def __init__(self, program_id, observing_time_fraction, field_ids,
                 filter_ids, internight_gap, n_visits_per_night,
                 intranight_gap, intranight_half_width,
                 nightly_priority='oldest', filter_choice='rotate'):

        self.program_id = program_id
        self.observing_time_fraction = observing_time_fraction
        self.field_ids = field_ids
        self.filter_ids = filter_ids

        self.internight_gap = internight_gap
        self.n_visits_per_night = n_visits_per_night
        self.intranight_gap = intranight_gap
        self.intranight_half_width = intranight_half_width

        self.nightly_priority = nightly_priority
        self.filter_choice = filter_choice

    def assign_nightly_requests(self, time, fields, **kwargs):

        # need a way to make this flexible without coding a new class for
        # every single way we could pick which fields to observe
        # some kind of config files?
        # or just instantiate with nightly_piority="", with
        # analagous functions to cadence.py?
        # selection functions:
        # nightly_priority='oldest":
        #    Prioritize fields for which it has been the longest
        #    time since the revisit
        # nightly_priority='radist':
        #    Prioritize fields closest to ra0
        # nightly_priority='decdist':
        #    Prioritize fields closest to dec0.  Give list to rotate???
        # nightly_priority='random'
        #    Choose fields at random

        # e.g.: fields assigned % n nights, rotate
        #               -> could handle this with cadence pars directly!
        #               but sensitive to state--get random initial set
        #               of fields  baked in
        #       observe oldest fields last observed
        #       observe fields on a fixed rotation cadence
        #       observe all fields
        #       observe fields nearest to an RA or Dec line
        #       observe randomly chosen fields

        # filters are given in filter_ids:
        # either a set of filters, or a fixed sequence
        # filter_choice = 'rotate':
        #   use one filter set per night, keyed by mjd % n_filters
        # filter_choice = 'sequence':
        #   use hard-coded sequence given in filter_ids

        # compute nightly altaz blocks and observability windows
        fields.compute_blocks(time)
        fields.compute_observability()

        n_filters = len(set(self.filter_ids))
        if self.filter_choice == 'rotate':
            night_index = np.floor(time.mjd % n_filters).astype(np.int)
            filter_ids_tonight = self.filter_ids[night_index]
        else:
            filter_ids_tonight = list(set(self.filter_ids))

        # Choose which fields will be observed

        obs_field_ids = fields.select_field_ids(
            # subtract an extra day since we are at the start of the night
            last_observed_range=[Time('2001-01-01').mjd,
                                 (time - self.internight_gap - 1 * u.day).mjd],
            program_id=self.program_id, filter_id=filter_ids_tonight,
            reducefunc=[np.min, np.min],  # we want oldest possible fields
            observable_hours_range=[1, 100.])

        # now form the intersection of observable fields and the OP fields
        pool_ids = obs_field_ids.intersection(self.field_ids)

        # how many fields can we observe tonight?
        n_requests = self.number_of_allowed_requests(time)

        n_fields = np.round(n_requests / self.n_visits_per_night)

        request_fields = fields.fields.loc[pool_ids]

        # sort request sets by chosen priority metric

        if self.nightly_priority == 'oldest':
            # now grab the top n_fields sorted by last observed date
            last_obs_keys = ['last_observed_{}_{}'.format(self.program_id, k)
                             for k in np.atleast_1d(filter_ids_tonight)]
            # make a new joint column
            oldest_obs = request_fields[last_obs_keys].apply(np.min, axis=1)
            oldest_obs.name = 'oldest_obs'
            request_fields = request_fields.join(oldest_obs)
            request_fields = request_fields.sort_values(
                by='oldest_obs').iloc[:n_fields]
            # last_obs_key = 'last_observed_{}_{}'.format(self.program_id,
            #                                            FILTER_NAME_TO_ID['r'])
            # request_fields = request_fields.sort_values(
            #    by=last_obs_key).iloc[:n_fields]
        elif self.nightly_priority == 'radist':
            assert ('ra0' in kwargs)
            # TODO: do I want spherical distance rather than ra distance?
            raise NotImplementedError
        elif self.nightly_priority == 'decdist':
            assert ('dec0' in kwargs)
            raise NotImplementedError
        elif self.nightly_priority == 'random':
            request_fields = request_fields.reindex(
                np.random.permutation(request_fields.index))[:n_fields]
        else:
            raise ValueError('requested prioritization scheme not found')

        # construct request sets: list of inputs to RequestPool.add_requests
        # scalar everything except field_ids

        if self.filter_choice == 'rotate':
            filter_sequence = [filter_ids_tonight for i in
                               range(self.n_visits_per_night)]
        elif self.filter_choice == 'sequence':
            assert(len(self.filter_ids) == self.n_visits_per_night)
            filter_sequence = self.filter_ids

        request_set = []
        # first visit
        request_set.append(
            {'program_id': self.program_id,
             'field_ids': request_fields.index,
             'filter_id': filter_sequence[0],
             'cadence_func': 'time_since_obs',
             'cadence_pars': {'ref_obs': 'last_observed',
                              'window_start': self.internight_gap.to(u.day).value,
                              # use a very large value here: gets added to
                              # last_obs.  remember that we reset each night
                              # anyway
                              'window_stop': (100 * u.year).to(u.day).value,
                              # TODO: do I want to specify this in some cases?
                              'prev_filter': 'any'},
             'priority': 1})
        # additional visits
        for i in range(self.n_visits_per_night - 1):
            request_set.append(
                {'program_id': self.program_id,
                 'field_ids': request_fields.index,
                 'filter_id': filter_sequence[i + 1],
                 'cadence_func': 'time_since_obs',
                 'cadence_pars': {'ref_obs': 'first_obs_tonight',
                                  'window_start': (i + 1) * (self.intranight_gap).to(u.day).value - self.intranight_half_width.to(u.day).value,
                                  'window_stop': (i + 1) * (self.intranight_gap).to(u.day).value + self.intranight_half_width.to(u.day).value,
                                  'prev_filter': 'any'},
                 'priority': 1})

        return request_set

    def number_of_allowed_requests(self, time):
        """ Count the (maximal) number of requests allowed for this program tonight."""

        # TODO: implement balancing of program observing time

        obs_time = approx_hours_of_darkness(
            time) * self.observing_time_fraction

        n_requests = (obs_time.to(u.min) /
                      (EXPOSURE_TIME + READOUT_TIME).to(u.min)).value
        return np.round(n_requests).astype(np.int)


# TODO: make these more functional, rather than hard-coding the behavior
class CollaborationObservingProgram(ObservingProgram):

    def __init__(self, field_ids,
                 filter_ids=[1, 1, 2, 2],
                 internight_gap=1. * u.day,
                 n_visits_per_night=4,
                 intranight_gap=60 * u.min,
                 intranight_half_width=20 * u.min,
                 nightly_priority='oldest',
                 filter_choice='sequence'):
        super(CollaborationObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['collaboration'], 0.4, field_ids, filter_ids,
            internight_gap, n_visits_per_night, intranight_gap,
            intranight_half_width, nightly_priority=nightly_priority,
            filter_choice=filter_choice)


class MSIPObservingProgram(ObservingProgram):

    def __init__(self, field_ids, filter_ids=[1, 2],
                 internight_gap=3. * u.day,
                 n_visits_per_night=2,
                 intranight_gap=60 * u.min,
                 intranight_half_width=20 * u.min,
                 nightly_priority='oldest',
                 filter_choice='rotate'):
        super(MSIPObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['MSIP'], 0.4, field_ids, filter_ids,
            internight_gap, n_visits_per_night, intranight_gap,
            intranight_half_width, nightly_priority=nightly_priority,
            filter_choice=filter_choice)


class CaltechObservingProgram(ObservingProgram):

    def __init__(self, field_ids,
                 filter_ids=[1, 2],
                 internight_gap=3. * u.day,
                 n_visits_per_night=3,
                 intranight_gap=90 * u.min,
                 intranight_half_width=20 * u.min,
                 nightly_priority='random',
                 filter_choice='rotate'):
        super(CaltechObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['Caltech'], 0.2, field_ids, filter_ids,
            internight_gap, n_visits_per_night, intranight_gap,
            intranight_half_width, nightly_priority=nightly_priority,
            filter_choice=filter_choice)
