import numpy as np
import astropy.units as u
import pytest
from ztf_sim.constants import slew_time, SETTLE_TIME


class TestSlewTime:

    def test_zero_angle_returns_zero(self):
        """Zero slew should take zero time (no settle time added)."""
        t = slew_time('ha', np.array([0.]) * u.deg)
        assert t[0] == 0.0 * u.second

    def test_small_angle_triangular_profile(self):
        """Very small slew uses triangular velocity profile (short accel/decel)."""
        t_small = slew_time('ha', np.array([0.1]) * u.deg)
        t_large = slew_time('ha', np.array([10.]) * u.deg)
        # Triangular profile: time grows as sqrt(angle); large slew is
        # proportionally faster per degree (trapezoidal)
        assert t_small[0] < t_large[0]

    def test_large_angle_trapezoidal_profile(self):
        """Large slew reaches vmax: time grows roughly linearly with angle."""
        t1 = slew_time('ha', np.array([30.]) * u.deg)
        t2 = slew_time('ha', np.array([60.]) * u.deg)
        # Doubling the angle should roughly double the time in the trapezoidal regime
        ratio = t2[0] / t1[0]
        assert 1.5 < ratio.value < 2.5

    @pytest.mark.parametrize('axis', ['ha', 'dec', 'dome'])
    def test_all_axes_return_quantity(self, axis):
        """Each axis returns an astropy Quantity in seconds."""
        t = slew_time(axis, np.array([5.]) * u.deg)
        assert t.unit.is_equivalent(u.second)

    def test_array_input(self):
        """Array of angles returns array of times with matching shape."""
        angles = np.array([0., 1., 5., 30.]) * u.deg
        t = slew_time('dec', angles)
        assert t.shape == angles.shape

    def test_nonzero_slew_adds_settle_time(self):
        """Nonzero slew includes SETTLE_TIME."""
        t = slew_time('ha', np.array([5.]) * u.deg)
        assert t[0] >= SETTLE_TIME
