import numpy as np
import pytest
from ztf_sim.magnitudes import limiting_mag


EXPOSURE_TIME = 30.  # seconds
SEEING = 2.0         # arcsec
SKY = 20.0           # mag/arcsec^2


class TestLimitingMag:

    def test_zenith_deeper_than_high_airmass(self):
        """Pointing at zenith (alt=90) gives deeper limit than low altitude."""
        lim_zenith = limiting_mag(EXPOSURE_TIME, SEEING, SKY, altitude=90.)
        lim_low = limiting_mag(EXPOSURE_TIME, SEEING, SKY, altitude=30.)
        assert lim_zenith > lim_low

    def test_longer_exposure_deeper(self):
        """Doubling exposure time increases limiting mag."""
        lim_short = limiting_mag(30., SEEING, SKY)
        lim_long = limiting_mag(120., SEEING, SKY)
        assert lim_long > lim_short

    def test_higher_snr_threshold_shallower(self):
        """Requiring higher SNR gives a shallower limiting magnitude."""
        lim_5sig = limiting_mag(EXPOSURE_TIME, SEEING, SKY, SNR=5.)
        lim_10sig = limiting_mag(EXPOSURE_TIME, SEEING, SKY, SNR=10.)
        assert lim_5sig > lim_10sig

    @pytest.mark.parametrize('filter_id', [1, 2, 3])
    def test_all_filters_return_finite(self, filter_id):
        """All three filters produce a finite limiting magnitude."""
        lim = limiting_mag(EXPOSURE_TIME, SEEING, SKY, filter_id=filter_id)
        assert np.isfinite(lim)

    def test_returns_positive_magnitude(self):
        """Result should be a positive AB magnitude."""
        lim = limiting_mag(EXPOSURE_TIME, SEEING, SKY)
        assert lim > 0
