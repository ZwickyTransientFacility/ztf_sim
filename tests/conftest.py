import pathlib
import pytest
from astropy.time import Time

REPO_ROOT = pathlib.Path(__file__).parent.parent
SIMS_DIR = REPO_ROOT / 'sims'


@pytest.fixture(scope='session')
def fields():
    from ztf_sim.Fields import Fields
    return Fields()


@pytest.fixture
def queue_config_path():
    return SIMS_DIR / 'survey_180501.json'


@pytest.fixture
def scheduler_config_path():
    return SIMS_DIR / 'example_scheduler_config.json'


@pytest.fixture
def obs_logger(tmp_path):
    from ztf_sim.ObsLogger import ObsLogger
    return ObsLogger('test_log', survey_start_time=Time('2018-01-01'),
                     output_path=str(tmp_path), clobber=True)


@pytest.fixture
def mock_obs_log():
    from unittest.mock import MagicMock
    import pandas as pd
    mock = MagicMock()
    # By default return empty DataFrame (field never observed)
    mock.select_last_observed_time_by_field.return_value = pd.DataFrame()
    return mock
