import numpy as np
import pandas as pd
import pytest
from ztf_sim.Fields import Fields


class TestFieldsLoading:

    def test_loads_without_error(self, fields):
        assert fields.fields is not None

    def test_field_count_reasonable(self, fields):
        """ZTF grid has ~1700 fields above dec=-32."""
        assert len(fields.fields) > 1000

    def test_index_is_field_id(self, fields):
        assert fields.fields.index.name == 'field_id'

    def test_columns_present(self, fields):
        for col in ['ra', 'dec', 'l', 'b', 'ecliptic_lon', 'ecliptic_lat', 'grid_id']:
            assert col in fields.fields.columns

    def test_dec_cut_applied(self, fields):
        """Fields below dec=-32 should have been dropped."""
        assert (fields.fields['dec'] >= -32).all()

    def test_grid_ids_valid(self, fields):
        """All grid_ids should be 0, 1, 2, or 3."""
        assert set(fields.fields['grid_id'].unique()).issubset({0, 1, 2, 3})


class TestFieldsSelection:

    def test_select_fields_dec_range(self, fields):
        cuts = fields.select_fields(dec_range=[30, 40])
        selected = fields.fields[cuts]
        assert (selected['dec'] >= 30).all()
        assert (selected['dec'] <= 40).all()
        assert len(selected) > 0

    def test_select_fields_ra_range(self, fields):
        cuts = fields.select_fields(ra_range=[180, 270])
        selected = fields.fields[cuts]
        assert (selected['ra'] >= 180).all()
        assert (selected['ra'] <= 270).all()

    def test_select_fields_grid_id(self, fields):
        cuts = fields.select_fields(grid_id=0)
        selected = fields.fields[cuts]
        assert (selected['grid_id'] == 0).all()
        assert len(selected) > 0

    def test_select_fields_abs_b_range(self, fields):
        cuts = fields.select_fields(abs_b_range=[30, 90])
        selected = fields.fields[cuts]
        assert (np.abs(selected['b']) >= 30).all()

    def test_b_range_and_abs_b_range_mutually_exclusive(self, fields):
        with pytest.raises(AssertionError):
            fields.select_fields(b_range=[10, 30], abs_b_range=[10, 30])

    def test_select_field_ids_returns_index(self, fields):
        result = fields.select_field_ids(dec_range=[30, 40])
        assert hasattr(result, 'name')  # pandas Index
        assert len(result) > 0

    def test_no_cuts_returns_all(self, fields):
        cuts = fields.select_fields()
        assert cuts.all()
