from __future__ import absolute_import

import numpy as np
from scipy.interpolate import interp1d
from .constants import *


# TODO: need to add i-band
def interp_R20_airmass(filter_id=2):
    """Returns function to interpolate electrons/sec of a 20th mag source as a function of altitude."""
    R20_file = '../data/R20_absorbed_ZTF{}.txt'.format(
        FILTER_ID_TO_NAME[filter_id])
    data = np.loadtxt(R20_file)
    alt = data[:, 0]
    R20 = data[:, 1]
    return interp1d(alt, R20)

R20_interp_alt = {1: interp_R20_airmass(filter_id=1),
                  2: interp_R20_airmass(filter_id=2),
                  3: interp_R20_airmass(filter_id=3)}


def limiting_mag(exposure_time, seeing_fwhm, sky_brightness,
                 filter_id=2, altitude=90., SNR=5.):
    """Calculate limiting magnitude."""

    npix = n_pixels(seeing_fwhm)
    Rsky = sky_electrons_per_pixel(sky_brightness, filter_id=filter_id)

    # sky limited case:
    Rstar = np.sqrt(SNR**2. * npix * Rsky / exposure_time)
    # TODO: decide if treatment of non-sky limited case is needed

    R20 = Rstar20(filter_id=filter_id, altitude=altitude,
                  aperture_cut=True, absorb=True)
    return 20. - 2.5 * np.log10(Rstar / R20)


def Rstar20(filter_id=2, altitude=90.,
            aperture_cut=True, absorb=True):
    """Compute electrons per second for a 20th mag source."""

    # make these arrays so we can handle input dataframes
    filter_id = np.atleast_1d(filter_id)
    altitude = np.atleast_1d(altitude)
    if len(altitude) == 1:
        altitude = np.ones(len(filter_id)) * altitude
    assert (len(filter_id) == len(altitude))
    R20 = np.zeros(len(filter_id))

    if (not aperture_cut) and (not absorb):
        # unabsorbed, no aperture cut: for calculating Rsky
        w1 = filter_id == 1
        R20[w1] = 123.73  # electrons/sec for 20th mag source
        w2 = filter_id == 2
        R20[w2] = 77.98
        w3 = filter_id == 3
        R20[w3] = 46.41
        if np.sum(w1) + np.sum(w2) + np.sum(w3) != len(R20):
            raise NotImplementedError

    elif aperture_cut and absorb:
        w1 = filter_id == 1
        R20[w1] = R20_interp_alt[1](altitude[w1])
        w2 = filter_id == 2
        R20[w2] = R20_interp_alt[2](altitude[w2])
        w3 = filter_id == 3
        R20[w3] = R20_interp_alt[3](altitude[w3])
        if np.sum(w1) + np.sum(w2) + np.sum(w3) != len(R20):
            raise NotImplementedError
    else:
        raise NotImplementedError

    return R20


def AB_to_Rstar(source_mag, filter_id=2, altitude=90.,
                aperture_cut=True, absorb=True):
    """Convert AB mag to electrons per second using a lookup table of 
    electrons per second for a 20th mag source."""

    R20 = Rstar20(filter_id=filter_id, altitude=altitude,
                  aperture_cut=True, absorb=True)

    return R20 * 10**(0.4 * (20. - source_mag))


def n_pixels(seeing_fwhm):
    """Calculate number of pixels in aperture extraction region.

    Uses the SNR-maximizing extraction radius
        (1.346 FWHM)--see LSST SNR doc eq. 19

    seeing_fwhm:  float
        seeing in arcsec
        """

    npix_extract = np.pi * (0.673 * seeing_fwhm / PIXEL_SCALE)**2.

    # don't return fractional pixels
    npix_extract = np.atleast_1d(np.round(npix_extract))

    w = npix_extract < 1.
    npix_extract[w] == 1

    return npix_extract


def sky_electrons_per_pixel(mag_per_sq_arcsec, filter_id=2):
    # returns electrons per pixel per second
    # area of one pixel in arcsec^2.
    pixarea = PIXEL_SCALE**2.
    mag_per_pix = mag_per_sq_arcsec - 2.5 * np.log10(pixarea)
    # could store R20 = R(20) and do R(m) = R(20) * 10**(0.4 * (20-m))
    return AB_to_Rstar(mag_per_pix, filter_id=filter_id, altitude=90.,
                       aperture_cut=False, absorb=False)
