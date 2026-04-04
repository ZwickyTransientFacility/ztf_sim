import pytest
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time
from ztf_sim.TelescopeStateMachine import TelescopeStateMachine
from ztf_sim.constants import P48_loc, FILTER_CHANGE_TIME, EXPOSURE_TIME


START_TIME = Time('2018-05-01 04:00:00', scale='utc', location=P48_loc)


@pytest.fixture
def tsm():
    return TelescopeStateMachine(
        current_time=START_TIME,
        historical_observability_year=None,  # no weather, faster
    )


class TestTelescopeStateMachineInit:

    def test_initial_state_is_ready(self, tsm):
        assert tsm.state == 'ready'

    def test_current_state_dict_keys(self, tsm):
        d = tsm.current_state_dict()
        for key in ['current_time', 'current_ha', 'current_dec',
                    'current_domeaz', 'current_filter_id',
                    'current_zenith_seeing', 'filters']:
            assert key in d


class TestFilterChange:

    def test_filter_change_advances_time(self, tsm):
        t0 = tsm.current_time.copy()
        tsm.start_filter_change(target_filter_id=1)  # from default r=2 to g=1
        assert tsm.current_time > t0

    def test_filter_change_time_is_correct(self, tsm):
        t0 = tsm.current_time.copy()
        tsm.start_filter_change(target_filter_id=1)
        dt = (tsm.current_time - t0).to(u.second)
        assert dt.value == pytest.approx(FILTER_CHANGE_TIME.to(u.second).value, rel=1e-3)

    def test_same_filter_no_time_advance(self, tsm):
        t0 = tsm.current_time.copy()
        tsm.start_filter_change(target_filter_id=tsm.current_filter_id)
        assert tsm.current_time == t0


class TestExposure:

    def test_exposure_advances_time(self, tsm):
        t0 = tsm.current_time.copy()
        # Need a target_skycoord set for process_exposure to compute HA/Az
        tsm.target_skycoord = coord.SkyCoord(180. * u.deg, 33. * u.deg)
        tsm.start_exposing(exposure_time=EXPOSURE_TIME)
        assert tsm.current_time > t0

    def test_exposure_time_matches(self, tsm):
        t0 = tsm.current_time.copy()
        tsm.target_skycoord = coord.SkyCoord(180. * u.deg, 33. * u.deg)
        tsm.start_exposing(exposure_time=EXPOSURE_TIME)
        dt = (tsm.current_time - t0).to(u.second)
        assert dt.value == pytest.approx(EXPOSURE_TIME.to(u.second).value, rel=1e-2)


class TestSlewAllowed:

    def test_below_horizon_rejected(self, tsm):
        # dec=-80 is well below Palomar's horizon at any HA
        target = coord.SkyCoord(0. * u.deg, -80. * u.deg)
        assert tsm.slew_allowed(target) is False

    def test_declination_below_limit_rejected(self, tsm):
        target = coord.SkyCoord(0. * u.deg, -36. * u.deg)
        assert tsm.slew_allowed(target) is False
