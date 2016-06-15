from constants import *
from utils import approx_hours_of_darkness


class ObservingProgram(object):

    def __init__(self, program_id, observing_time_fraction, field_ids,
                 internight_gap, n_visits_per_night, intranight_gap_range):

        self.program_id = program_id
        self.observing_time_fraction = observing_time_fraction
        self.field_ids = field_ids

        self.internight_gap = internight_gap
        self.n_visits_per_night = n_visits_per_night
        self.intranight_gap_range = intranight_gap_range

    def assign_nightly_requests(self, time, fields):

        # need a way to make this flexible without coding a new class for
        # every single way we could pick which fields to observe
        # some kind of config files?  or just instantiate with how="", with
        # analagous functions to cadence.py?

        # e.g.: fields assigned % n nights, rotate
        #               -> could handle this with cadence pars directly!
        #               but sensitive to state--get random initial set
        #               of fields  baked in
        #       observe oldest fields last observed
        #       observe all fields

        # just subclass, for now
        return self._assign_nightly_requests(time, fields)

    def _assign_nightly_requests(self, time, fields):

        # two visits per night, separated by an hour.  Prioritize
        # fields for which it has been the longest time since the revisit

        # find fields in current sample at least N days since last obs
        # internight_gap = 3. * u.day

        # compute nightly altaz blocks and observability windows
        fields.compute_blocks(time)
        fields.compute_observability()

        obs_field_ids = fields.select_field_ids(last_observed_range=# subtract an extra day since
                                                # we are at the start of the
                                                # night
                                                [Time('2001-01-01').mjd,
                                                 (time - self.internight_gap - 1 * u.day).mjd],
                                                program_id=self.program_id,
                                                filter_id=FILTER_NAME_TO_ID[
                                                    'r'],
                                                observable_hours_range=[1, 100.])
        # now form the intersection of observable fields and the OP fields
        pool_ids = obs_field_ids.intersection(self.field_ids)

        # how many fields can we observe tonight?
        n_requests = self.number_of_allowed_requests(time)
        # n_visits_per_night = 2

        n_fields = np.round(n_requests / self.n_visits_per_night)

        request_fields = fields.fields.loc[pool_ids]

        # now grab the top n_fields sorted by last observed date
        last_obs_key = 'last_observed_{}_{}'.format(self.program_id,
                                                    FILTER_NAME_TO_ID['r'])
        request_fields = request_fields.sort_values(
            by=last_obs_key).iloc[:n_fields]

        # construct request sets: list of inputs to RequestPool.add_requests
        # scalar everything except field_ids

        request_set = []
        # first visit
        request_set.append(
            {'program_id': self.program_id,
             'field_ids': request_fields.index,
             'filter_id': FILTER_NAME_TO_ID['r'],
             'cadence_func': 'time_since_last_obs',
             'cadence_pars': {'window_start': self.internight_gap.to(u.day).value,
                              # inf causes problems
                              'window_stop': (100 * u.year).to(u.day).value,
                              'prev_filter': 'any'},
             'priority': 1})
        # second visit
        for i in range(self.n_visits_per_night - 1):
            request_set.append(
                {'program_id': self.program_id,
                 'field_ids': request_fields.index,
                 'filter_id': FILTER_NAME_TO_ID['r'],
                 'cadence_func': 'time_since_last_obs',
                 'cadence_pars': {'window_start': (self.intranight_gap_range[0]).to(u.day).value,
                                  'window_stop': (self.intranight_gap_range[1]).to(u.day).value,
                                  'prev_filter': 'any'},
                 'priority': 1})

        return request_set

    def number_of_allowed_requests(self, time):
        """ Count the (maximal) number of requests allowed for this program tonight."""

        obs_time = approx_hours_of_darkness(
            time) * self.observing_time_fraction

        n_requests = (obs_time.to(u.min) /
                      (EXPOSURE_TIME + READOUT_TIME).to(u.min)).value
        return np.round(n_requests).astype(np.int)


# TODO: make these more functional, rather than hard-coding the behavior
class CollaborationObservingProgram(ObservingProgram):

    def __init__(self, field_ids,
                 internight_gap=1. * u.day,
                 n_visits_per_night=10,
                 intranight_gap_range=[
                     40 * u.min, 100 * u.min]):
        super(CollaborationObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['collaboration'], 0.4, field_ids,
            internight_gap, n_visits_per_night, intranight_gap_range)


class MSIPObservingProgram(ObservingProgram):

    def __init__(self, field_ids,
                 internight_gap=3. * u.day,
                 n_visits_per_night=2,
                 intranight_gap_range=[
                     40 * u.min, 100 * u.min]):
        super(MSIPObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['MSIP'], 0.4, field_ids,
            internight_gap, n_visits_per_night, intranight_gap_range)


class CaltechObservingProgram(ObservingProgram):

    def __init__(self, field_ids,
                 internight_gap=3. * u.day,
                 n_visits_per_night=2,
                 intranight_gap_range=[
                     40 * u.min, 100 * u.min]):
        super(CaltechObservingProgram, self).__init__(
            PROGRAM_NAME_TO_ID['Caltech'], 0.2, field_ids)

    def _assign_nightly_requests(self, time, fields):
        raise NotImplementedError
