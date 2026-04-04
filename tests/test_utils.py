import numpy as np
import pytest
import astropy.units as u
from astropy.time import Time
from ztf_sim.utils import (
    block_index, nightly_blocks, altitude_to_airmass,
    RA_to_HA, HA_to_RA, airmass_to_altitude,
)
from ztf_sim.constants import P48_loc


TEST_TIME = Time('2018-05-01 04:15:00', scale='utc', location=P48_loc)


class TestBlockIndex:

    def test_returns_integer_array(self):
        result = block_index(TEST_TIME)
        assert result.dtype in (np.int32, np.int64)

    def test_increases_with_time(self):
        t1 = TEST_TIME
        t2 = TEST_TIME + 30 * u.minute
        b1 = block_index(t1)
        b2 = block_index(t2)
        assert b2[0] >= b1[0]

    def test_same_block_within_1min(self):
        t1 = TEST_TIME
        t2 = TEST_TIME + 1 * u.minute
        b1 = block_index(t1)
        b2 = block_index(t2)
        assert b1[0] == b2[0]


class TestNightlyBlocks:

    def test_returns_two_arrays(self):
        blocks, times = nightly_blocks(TEST_TIME)
        assert len(blocks) > 0
        assert len(times) == len(blocks)

    def test_block_count_reasonable(self):
        """A May night at Palomar is roughly 8-9 hours -> 16-18 30-min blocks."""
        blocks, _ = nightly_blocks(TEST_TIME)
        assert 10 < len(blocks) < 30

    def test_times_are_monotonically_increasing(self):
        _, times = nightly_blocks(TEST_TIME)
        mjds = times.mjd
        assert (np.diff(mjds) > 0).all()


class TestAirmass:

    def test_zenith_airmass_is_one(self):
        result = altitude_to_airmass(90.)
        assert result == pytest.approx(1.0, rel=1e-3)

    def test_airmass_increases_toward_horizon(self):
        x_high = altitude_to_airmass(60.)
        x_low = altitude_to_airmass(30.)
        assert x_low > x_high

    def test_airmass_altitude_roundtrip(self):
        """altitude -> airmass -> altitude should round-trip."""
        alt_in = 45.
        airmass = altitude_to_airmass(alt_in)
        alt_out = airmass_to_altitude(airmass).value
        assert alt_out == pytest.approx(alt_in, rel=1e-3)


class TestRAHA:

    def test_ra_to_ha_and_back(self):
        """HA_to_RA(RA_to_HA(ra)) should recover ra."""
        ra_in = 180. * u.deg
        ha = RA_to_HA(ra_in, TEST_TIME)
        ra_out = HA_to_RA(ha, TEST_TIME)
        assert ra_out.value == pytest.approx(ra_in.value, abs=0.01)

    def test_ha_to_ra_and_back(self):
        """RA_to_HA(HA_to_RA(ha)) should recover ha."""
        ha_in = 45. * u.deg
        ra = HA_to_RA(ha_in, TEST_TIME)
        ha_out = RA_to_HA(ra, TEST_TIME)
        assert ha_out.value == pytest.approx(ha_in.value, abs=0.01)
