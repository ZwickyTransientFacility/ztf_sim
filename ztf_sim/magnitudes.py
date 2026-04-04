"""Utilities for magnitude conversions"""

import numpy as np
from scipy.interpolate import interp1d
from .constants import BASE_DIR, FILTER_ID_TO_NAME, PIXEL_SCALE




def interp_R20_airmass(filter_id=2):
    """Build an interpolator for the 20th-magnitude electron rate vs. airmass.

    Reads a pre-computed lookup table from
    ``data/R20_absorbed_ZTF{filter}.txt`` and returns a callable that maps
    pointing altitude (degrees) to the electron rate (e⁻/s) for a 20th AB
    magnitude point source observed through the full ZTF optical path.

    Parameters
    ----------
    filter_id : int, optional
        Filter identifier: 1 = g, 2 = r, 3 = i. Default is 2.

    Returns
    -------
    scipy.interpolate.interp1d
        Callable that takes altitude in degrees and returns electron rate
        (e⁻/s) for a 20th AB mag source.
    """
    R20_file = BASE_DIR + '../data/R20_absorbed_ZTF{}.txt'.format(
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
    """Compute the point-source limiting AB magnitude for a given exposure.

    Assumes the sky-limited regime. Uses the SNR-maximising extraction
    aperture (radius = 1.346 × FWHM, see `n_pixels`).

    Parameters
    ----------
    exposure_time : float
        Exposure duration in seconds.
    seeing_fwhm : float or array-like
        Seeing FWHM at the pointing altitude in arcseconds.
    sky_brightness : float or array-like
        Sky surface brightness in mag arcsec⁻².
    filter_id : int or array-like, optional
        Filter identifier(s): 1 = g, 2 = r, 3 = i. Default is 2.
    altitude : float or array-like, optional
        Pointing altitude in degrees. Default is 90 (zenith).
    SNR : float, optional
        Required signal-to-noise ratio. Default is 5.

    Returns
    -------
    float or numpy.ndarray
        5σ limiting AB magnitude.
    """

    npix = n_pixels(seeing_fwhm)
    Rsky = sky_electrons_per_pixel(sky_brightness, filter_id=filter_id)

    # sky limited case:
    Rstar = np.sqrt(SNR**2. * npix * Rsky / exposure_time)

    R20 = Rstar20(filter_id=filter_id, altitude=altitude,
                  aperture_cut=True, absorb=True)
    return 20. - 2.5 * np.log10(Rstar / R20)


def Rstar20(filter_id=2, altitude=90.,
            aperture_cut=True, absorb=True):
    """Compute the electron rate for a 20th AB magnitude point source.

    Parameters
    ----------
    filter_id : int or array-like
        Filter identifier(s): 1 = g, 2 = r, 3 = i. Default is 2.
    altitude : float or array-like, optional
        Pointing altitude in degrees. Scalar or same length as *filter_id*.
        Default is 90 (zenith).
    aperture_cut : bool, optional
        If ``True``, apply the finite-aperture correction (default). Must be
        paired with ``absorb=True``; the ``True``/``False`` and
        ``False``/``True`` combinations raise ``NotImplementedError``.
    absorb : bool, optional
        If ``True``, include atmospheric absorption (default). See note above.

    Returns
    -------
    numpy.ndarray
        Electron rate (e⁻/s) with the same length as *filter_id*.

    Raises
    ------
    NotImplementedError
        If an unknown filter ID is encountered.

    Notes
    -----
    ``aperture_cut=True, absorb=True``: reads altitude-dependent values from
    ``data/R20_absorbed_ZTF{filter}.txt`` via ``interp_R20_airmass``.

    ``aperture_cut=False, absorb=False``: returns fixed zenith values
    (g = 123.73, r = 77.98, i = 46.41 e⁻/s). Used for sky electron
    calculations.
    """

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
    """Convert an AB magnitude to an electron count rate.

    Parameters
    ----------
    source_mag : float or array-like
        AB magnitude of the source.
    filter_id : int or array-like, optional
        Filter identifier(s): 1 = g, 2 = r, 3 = i. Default is 2.
    altitude : float or array-like, optional
        Pointing altitude in degrees. Default is 90 (zenith).
    aperture_cut : bool, optional
        Passed to `Rstar20`. Default is ``True``.
    absorb : bool, optional
        Passed to `Rstar20`. Default is ``True``.

    Returns
    -------
    numpy.ndarray
        Electron rate (e⁻/s) for the source.

    Notes
    -----
    Scaling relation::

        Rstar = R20 * 10^(0.4 * (20 - source_mag))
    """

    R20 = Rstar20(filter_id=filter_id, altitude=altitude,
                  aperture_cut=True, absorb=True)

    return R20 * 10**(0.4 * (20. - source_mag))


def n_pixels(seeing_fwhm):
    """Count the pixels within the SNR-maximising extraction aperture.

    Parameters
    ----------
    seeing_fwhm : float or array-like
        Seeing FWHM at the pointing altitude in arcseconds.

    Returns
    -------
    numpy.ndarray
        Number of pixels (rounded, minimum 1) within the extraction aperture.
    """

    npix_extract = np.pi * (0.673 * seeing_fwhm / PIXEL_SCALE)**2.

    # don't return fractional pixels
    npix_extract = np.atleast_1d(np.round(npix_extract))

    w = npix_extract < 1.
    npix_extract[w] == 1

    return npix_extract


def sky_electrons_per_pixel(mag_per_sq_arcsec, filter_id=2):
    """Convert sky surface brightness to electron rate per pixel per second.

    Parameters
    ----------
    mag_per_sq_arcsec : float or array-like
        Sky surface brightness in AB mag arcsec⁻².
    filter_id : int or array-like, optional
        Filter identifier(s): 1 = g, 2 = r, 3 = i. Default is 2.

    Returns
    -------
    float or numpy.ndarray
        Sky background electron rate in e⁻ pixel⁻¹ s⁻¹.
    """
    # area of one pixel in arcsec^2.
    pixarea = PIXEL_SCALE**2.
    mag_per_pix = mag_per_sq_arcsec - 2.5 * np.log10(pixarea)
    # could store R20 = R(20) and do R(m) = R(20) * 10**(0.4 * (20-m))
    return AB_to_Rstar(mag_per_pix, filter_id=filter_id, altitude=90.,
                       aperture_cut=False, absorb=False)
