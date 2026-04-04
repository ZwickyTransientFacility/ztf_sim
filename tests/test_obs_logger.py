import os
import pytest
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time
from ztf_sim.ObsLogger import ObsLogger
from ztf_sim.constants import P48_loc


START_TIME = Time('2018-05-01', scale='utc')
OBS_TIME = Time('2018-05-01 04:30:00', scale='utc', location=P48_loc)


def _make_state():
    return {
        'current_time': OBS_TIME + 30. * u.second,
        'current_zenith_seeing': 2.0 * u.arcsec,
    }


def _make_request():
    return {
        'request_id': 1,
        'target_program_id': 1,
        'target_field_id': 635,
        'target_ra': 180.0,
        'target_dec': 30.0,
        'target_filter_id': 2,
        'target_exposure_time': 30. * u.second,
        'target_sky_brightness': 20.0,
        'target_limiting_mag': 20.5,
        'target_total_requests_tonight': 5,
        'target_metric_value': 1.0,
        'target_subprogram_name': 'test',
    }


class TestObsLoggerSetup:

    def test_creates_db_file(self, tmp_path):
        logger = ObsLogger('mylog', survey_start_time=START_TIME,
                           output_path=str(tmp_path), clobber=True)
        assert (tmp_path / 'mylog.db').exists()

    def test_history_is_empty_on_init(self, obs_logger):
        assert len(obs_logger.history) == 0


class TestLogPointing:

    def test_log_pointing_adds_row(self, obs_logger):
        obs_logger.log_pointing(_make_state(), _make_request())
        assert len(obs_logger.history) == 1

    def test_log_pointing_records_field_id(self, obs_logger):
        obs_logger.log_pointing(_make_state(), _make_request())
        assert obs_logger.history.iloc[0]['fieldID'] == 635

    def test_log_pointing_records_filter(self, obs_logger):
        obs_logger.log_pointing(_make_state(), _make_request())
        assert obs_logger.history.iloc[0]['filter'] == 'r'

    def test_select_last_observed_returns_mjd(self, obs_logger):
        obs_logger.log_pointing(_make_state(), _make_request())
        result = obs_logger.select_last_observed_time_by_field(
            field_ids={635}, program_ids=[1], subprogram_names=['test'])
        assert 635 in result.index
        assert np.isfinite(result.loc[635, 'expMJD'])

    def test_two_pointings_accumulate(self, obs_logger):
        obs_logger.log_pointing(_make_state(), _make_request())
        obs_logger.log_pointing(_make_state(), _make_request())
        assert len(obs_logger.history) == 2
