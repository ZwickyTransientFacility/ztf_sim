import numpy as np
import pandas as pd
import pytest
import astropy.units as u
from astropy.time import Time
from unittest.mock import MagicMock
from ztf_sim.cadence import no_cadence, enough_gap_since_last_obs


NOW_MJD = 58200.5  # arbitrary reference MJD
NOW = Time(NOW_MJD, format='mjd')
GAP_MIN = 30.0  # required gap in minutes
GAP_DAYS = GAP_MIN / 1440.


def _make_df(field_ids, gap_min=GAP_MIN):
    return pd.DataFrame({
        'field_id': field_ids,
        'program_id': [1] * len(field_ids),
        'subprogram_name': ['test'] * len(field_ids),
        'intranight_gap_min': [gap_min] * len(field_ids),
    })


def _make_state():
    return {'current_time': NOW}


def _mock_obs_log_with(field_id, last_obs_mjd):
    """Return a mock obs_log that reports field_id was last observed at last_obs_mjd."""
    mock = MagicMock()
    mock.select_last_observed_time_by_field.return_value = pd.DataFrame(
        {'expMJD': [last_obs_mjd]}, index=pd.Index([field_id], name='fieldID'))
    return mock


def _mock_obs_log_empty():
    mock = MagicMock()
    mock.select_last_observed_time_by_field.return_value = pd.DataFrame()
    return mock


class TestNoCadence:
    def test_always_returns_true(self):
        assert no_cadence() is True

    def test_ignores_arguments(self):
        assert no_cadence(1, 2, 'anything') is True


class TestEnoughGapSinceLastObs:

    def test_never_observed_field_always_passes(self):
        df = _make_df([635])
        result = enough_gap_since_last_obs(df, _make_state(), _mock_obs_log_empty())
        assert result.all()

    def test_recently_observed_field_fails(self):
        """Field observed 10 min ago fails a 30-min gap requirement."""
        df = _make_df([635])
        recent_mjd = NOW_MJD - 10. / 1440.
        result = enough_gap_since_last_obs(
            df, _make_state(), _mock_obs_log_with(635, recent_mjd))
        assert not result.any()

    def test_old_enough_observation_passes(self):
        """Field observed 60 min ago passes a 30-min gap requirement."""
        df = _make_df([635])
        old_mjd = NOW_MJD - 60. / 1440.
        result = enough_gap_since_last_obs(
            df, _make_state(), _mock_obs_log_with(635, old_mjd))
        assert result.all()

    def test_exactly_at_gap_passes(self):
        """Field observed exactly GAP_MIN ago is eligible (>=)."""
        df = _make_df([635])
        exact_mjd = NOW_MJD - GAP_DAYS
        result = enough_gap_since_last_obs(
            df, _make_state(), _mock_obs_log_with(635, exact_mjd))
        assert result.all()

    def test_mixed_dataframe(self):
        """One passing and one failing field returns correct per-row mask."""
        df = _make_df([635, 636])
        # field 635: observed 60 min ago (passes); field 636: 5 min ago (fails)
        def side_effect(field_ids, program_ids, subprogram_names):
            ids = list(field_ids)
            rows = []
            for fid in ids:
                if fid == 635:
                    rows.append(NOW_MJD - 60. / 1440.)
                else:
                    rows.append(NOW_MJD - 5. / 1440.)
            return pd.DataFrame({'expMJD': rows},
                                index=pd.Index(ids, name='fieldID'))
        mock = MagicMock()
        mock.select_last_observed_time_by_field.side_effect = side_effect
        result = enough_gap_since_last_obs(df, _make_state(), mock)
        assert bool(result.iloc[0]) is True
        assert bool(result.iloc[1]) is False
