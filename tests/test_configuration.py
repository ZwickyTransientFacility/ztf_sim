import json
import pathlib
import pytest
from ztf_sim.configuration import SchedulerConfiguration, QueueConfiguration

SIMS_DIR = pathlib.Path(__file__).parent.parent / 'sims'


class TestSchedulerConfiguration:

    def test_valid_config_loads(self):
        cfg = SchedulerConfiguration(SIMS_DIR / 'example_scheduler_config.json')
        assert 'queues' in cfg.config

    def test_missing_queues_key_raises(self, tmp_path):
        bad = tmp_path / 'bad.json'
        bad.write_text(json.dumps({'run_name': 'test'}))
        with pytest.raises(ValueError, match='queues'):
            SchedulerConfiguration(bad)

    def test_missing_default_queue_raises(self, tmp_path):
        bad = tmp_path / 'bad.json'
        bad.write_text(json.dumps({
            'queues': [{'queue_name': 'not_default', 'config_file': 'x.json'}]
        }))
        with pytest.raises(ValueError, match='default'):
            SchedulerConfiguration(bad)

    def test_build_queue_configs_returns_dict(self):
        cfg = SchedulerConfiguration(SIMS_DIR / 'example_scheduler_config.json')
        queue_configs = cfg.build_queue_configs()
        assert 'default' in queue_configs


class TestQueueConfiguration:

    def test_valid_config_loads(self):
        cfg = QueueConfiguration(SIMS_DIR / 'survey_180501.json')
        assert 'observing_programs' in cfg.config

    def test_fractions_not_summing_to_one_raises(self, tmp_path):
        bad = tmp_path / 'bad_fracs.json'
        bad.write_text(json.dumps({
            'queue_manager': 'greedy',
            'observing_programs': [
                {
                    'program_name': 'MSIP',
                    'subprogram_name': 'test',
                    'program_pi': 'pi',
                    'program_observing_fraction': 0.5,
                    'subprogram_fraction': 0.5,  # 0.5*0.5 = 0.25, not 1.0
                    'active_months': 'all',
                    'filter_ids': [2],
                    'internight_gap_days': 3,
                    'n_visits_per_night': 1,
                    'filter_choice': 'rotate',
                    'field_ids': [635],
                }
            ]
        }))
        with pytest.raises(ValueError, match='sum'):
            QueueConfiguration(bad)

    def test_unknown_program_name_raises(self, tmp_path):
        bad = tmp_path / 'bad_prog.json'
        bad.write_text(json.dumps({
            'queue_manager': 'greedy',
            'observing_programs': [
                {
                    'program_name': 'NONEXISTENT_PROGRAM',
                    'subprogram_name': 'test',
                    'program_pi': 'pi',
                    'program_observing_fraction': 1.0,
                    'subprogram_fraction': 1.0,
                    'active_months': 'all',
                    'filter_ids': [2],
                    'internight_gap_days': 3,
                    'n_visits_per_night': 1,
                    'filter_choice': 'rotate',
                    'field_ids': [635],
                }
            ]
        }))
        with pytest.raises(ValueError):
            QueueConfiguration(bad)

    def test_list_queue_manager_skips_fraction_check(self, tmp_path):
        """list queue_manager type does not require fractions to sum to 1."""
        cfg_data = tmp_path / 'list_queue.json'
        cfg_data.write_text(json.dumps({
            'queue_manager': 'list',
            'observing_programs': []
        }))
        cfg = QueueConfiguration(cfg_data)
        assert cfg.config['queue_manager'] == 'list'
