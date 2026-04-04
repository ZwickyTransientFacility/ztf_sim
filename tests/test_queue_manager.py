import pathlib
import pytest

gurobipy = pytest.importorskip('gurobipy')

from ztf_sim.QueueManager import ListQueueManager, GreedyQueueManager, GurobiQueueManager
from ztf_sim.configuration import QueueConfiguration

SIMS_DIR = pathlib.Path(__file__).parent.parent / 'sims'


@pytest.fixture
def list_queue_config(tmp_path):
    import json
    cfg = tmp_path / 'list.json'
    cfg.write_text(json.dumps({
        'queue_manager': 'list',
        'observing_programs': [],
        'targets': [{'field_id': 635, 'program_id': 1,
                     'subprogram_name': 'test', 'filter_id': 2,
                     'program_pi': 'pi'}]
    }))
    return QueueConfiguration(cfg)


@pytest.fixture
def queue_config():
    return QueueConfiguration(SIMS_DIR / 'survey_180501.json')


class TestListQueueManager:

    def test_load_list_queue(self, list_queue_config):
        lqm = ListQueueManager('test', list_queue_config)
        lod = [
            {'field_id': 635, 'program_id': 1, 'subprogram_name': 'test',
             'filter_id': 2, 'program_pi': 'pi'},
            {'field_id': 680, 'program_id': 1, 'subprogram_name': 'test',
             'filter_id': 1, 'program_pi': 'pi'},
        ]
        lqm.load_list_queue(lod)
        assert len(lqm.queue) == 2

    def test_load_list_queue_append(self, list_queue_config):
        lqm = ListQueueManager('test', list_queue_config)
        lod = [{'field_id': 635, 'program_id': 1, 'subprogram_name': 'test',
                'filter_id': 2, 'program_pi': 'pi'}]
        lqm.load_list_queue(lod)
        lqm.load_list_queue(lod, append=True)
        assert len(lqm.queue) == 2

    def test_load_list_queue_replace(self, list_queue_config):
        lqm = ListQueueManager('test', list_queue_config)
        lod = [{'field_id': 635, 'program_id': 1, 'subprogram_name': 'test',
                'filter_id': 2, 'program_pi': 'pi'}]
        lqm.load_list_queue(lod)
        lqm.load_list_queue(lod)  # replace, not append
        assert len(lqm.queue) == 1


class TestGreedyQueueManager:

    def test_instantiates(self, queue_config):
        gqm = GreedyQueueManager('test', queue_config)
        assert gqm is not None


class TestGurobiQueueManager:

    def test_instantiates(self, queue_config):
        gqm = GurobiQueueManager('test', queue_config)
        assert gqm is not None
