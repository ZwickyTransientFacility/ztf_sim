import pytest

gurobipy = pytest.importorskip('gurobipy')

from ztf_sim.optimize import night_optimize, tsp_optimize


class TestOptimizeImports:

    def test_night_optimize_callable(self):
        assert callable(night_optimize)

    def test_tsp_optimize_callable(self):
        assert callable(tsp_optimize)
