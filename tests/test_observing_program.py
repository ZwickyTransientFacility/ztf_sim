import pytest
import astropy.units as u
from ztf_sim.ObservingProgram import ObservingProgram


def _make_op(**kwargs):
    defaults = dict(
        program_id=1,
        subprogram_name='test',
        program_pi='pi',
        program_observing_time_fraction=1.0,
        subprogram_fraction=1.0,
        field_ids=[635],
        filter_ids=[2],
        internight_gap=3 * u.day,
        intranight_gap=30 * u.minute,
        n_visits_per_night=1,
        field_selection_function=None,
        active_months='all',
        filter_choice='rotate',
    )
    defaults.update(kwargs)
    return ObservingProgram(**defaults)


class TestObservingProgramInit:

    def test_accepts_field_ids(self):
        op = _make_op(field_ids=[635, 680], field_selection_function=None)
        assert op.field_ids == [635, 680]

    def test_accepts_field_selection_function(self):
        op = _make_op(field_ids=None, field_selection_function='MSIPSurvey')
        assert op.field_selection_function == 'MSIPSurvey'

    def test_raises_if_both_provided(self):
        with pytest.raises(AssertionError):
            _make_op(field_ids=[635], field_selection_function='MSIPSurvey')

    def test_raises_if_neither_provided(self):
        with pytest.raises(AssertionError):
            _make_op(field_ids=None, field_selection_function=None)

    def test_active_months_all(self):
        op = _make_op(active_months='all')
        assert op.active_months == 'all'

    def test_active_months_list(self):
        op = _make_op(active_months=[4, 5, 6])
        assert 4 in op.active_months

    def test_filter_choice_stored(self):
        op = _make_op(filter_choice='sequence')
        assert op.filter_choice == 'sequence'

    def test_exposure_time_stored(self):
        op = _make_op(exposure_time=60 * u.second)
        assert op.exposure_time.to(u.second).value == pytest.approx(60.)
