"""
@author: yuhanyao
"""
from glob import glob
import os
import logging
import astropy.constants as const
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import get_sun
from ..Fields import Fields
from ..constants import BASE_DIR

logger = logging.getLogger(__name__)



def get_theta_given_phi(a, b, c, phi):
    """Solve for the polar angle θ on the plane ax + by + cz = 0 at azimuth φ.

    Given the plane equation ``a*x + b*y + c*z = 0`` (perpendicular to the
    Sun) and a fixed azimuthal angle φ, computes the polar angle θ such that
    the unit-sphere point ``(sin θ cos φ, sin θ sin φ, cos θ)`` lies on the
    plane.

    Parameters
    ----------
    a, b, c : float
        Coefficients of the plane equation (direction cosines of the Sun
        in Cartesian coordinates).
    phi : float
        Azimuthal angle in radians.

    Returns
    -------
    float
        Polar angle θ in radians, adjusted to lie in [0, π].
    """
    a_prime = a * np.cos(phi)
    b_prime = b * np.sin(phi)
    theta_ = np.arctan(-1 * c / (a_prime + b_prime))
    if theta_ < 0:
        theta = theta_+np.pi
    else:
        theta = theta_
    return theta


def cos_in_sphere(theta1, phi1, theta2, phi2):
    """Compute the cosine of the great-circle angular distance between two points.

    Uses the spherical law of cosines:
    ``cos(d) = cos(θ2)*cos(θ1) + sin(θ2)*sin(θ1)*cos(φ1 − φ2)``.

    Parameters
    ----------
    theta1, theta2 : float
        Polar (co-latitude) angles in radians.
    phi1, phi2 : float
        Azimuthal angles in radians.

    Returns
    -------
    float
        Cosine of the angular separation between the two points.
    """
    term1 = np.cos(theta2) * np.cos(theta1)
    term2 = np.sin(theta2) * np.sin(theta1) * np.cos(phi1-phi2)
    cos_rad = term1 + term2
    return cos_rad

def load_SRG_pointing(tnow):
    """Load the pre-computed SRG pointing plan for the current night.

    Reads all CSV files from ``ztf_survey_configuration/srg_plan/``,
    concatenates them, and filters to rows whose MJD falls within the
    current night (``[tnow.mjd, tnow.mjd + 1)``).

    Parameters
    ----------
    tnow : astropy.time.Time
        Any time within the night of interest.

    Returns
    -------
    alphas : numpy.ndarray of float
        Right ascensions of SRG pointing centres in degrees.
    deltas : numpy.ndarray of float
        Declinations of SRG pointing centres in degrees.

    Raises
    ------
    FileNotFoundError
        If no CSV files are found in the SRG plan directory.
    ValueError
        If no rows in the plan files cover the requested date.
    """
    csvfiles = glob(f"{BASE_DIR}../../ztf_survey_configuration/srg_plan/*.csv")
    if len(csvfiles) == 0:
        raise FileNotFoundError('No detailed SRG pointing plans found')
    csvfiles = np.array(csvfiles)
    csvfiles = csvfiles[np.argsort(csvfiles)]
    nfiles = len(csvfiles)
    for i in range(nfiles):
        if i==0:
            tb = pd.read_csv(csvfiles[i])
        else:
            tb = pd.concat([tb, pd.read_csv(csvfiles[i])])
    
    # TODO: passing in the desired date directly would be more efficient
    ix = (tb["MJD"]>=tnow.mjd)&(tb["MJD"]<(tnow.mjd+1))
    if np.sum(ix) == 0:
        raise ValueError(f'No detailed SRG pointing plans found for date {tnow.mjd}')
    tb = tb[ix]
    alphas = tb["RA"].values
    deltas = tb["Dec"].values
    return alphas, deltas

def SRG_pointing(tnow):
    """Compute the SRG/eROSITA scanning footprint for the current night.

    Models the SRG great-circle scan as a spiral perpendicular to the Sun,
    uniformly sampling ~360 points around the full circle. Near the spring
    or autumn equinox (when the scan passes through the poles), a finer
    step size is used to ensure adequate coverage.

    Parameters
    ----------
    tnow : astropy.time.Time
        Any time within the night of interest. The Sun position is evaluated
        at Palomar midnight (UTC 07:00).

    Returns
    -------
    alphas : numpy.ndarray of float
        Right ascensions of SRG pointing centres in degrees.
    decs : numpy.ndarray of float
        Declinations of SRG pointing centres in degrees.

    Notes
    -----
    SRG observing strategy: 1 sq deg FoV, ~1 sq deg per day drift, one full
    great-circle scan every ~4 hours. The scan plane is perpendicular to the
    Sun direction: ``a*x + b*y + c*z = 0`` where (a, b, c) are the Cartesian
    coordinates of the Sun.

    References
    ----------
    SRG real-time tracking: http://plan.srg.cosmos.ru/monthplan/tracking
    """
    # Palomar midnight = UTC 7am
    midnight_mjd = np.floor(tnow.mjd) + 7/24.
    sunpos = get_sun(Time(midnight_mjd, format='mjd'))
    sun_ra_rad = sunpos.ra.rad
    sub_dec_rad = sunpos.dec.rad
    
    sun_theta = np.pi/2 - sub_dec_rad
    sun_phi = sun_ra_rad
    a = np.sin(sun_theta)*np.cos(sun_phi)
    b = np.sin(sun_theta)*np.sin(sun_phi)
    c = np.cos(sun_theta)
    # the plane that is perpendicular to the Sun is 
    # ax + by + cz = 0
    phis = [0]
    thetas = [get_theta_given_phi(a, b, c, phis[0])]
    sep_degs = []

    # if we are near the spring or fall equinox (Sun RA = 360 or 180 deg), 
    # the plane passes through the poles and the default sampling is not enough
    is_spring_equinox = (sunpos.ra.degree > 357) or (sunpos.ra.degree < 3)
    is_fall_equinox = (sunpos.ra.degree > 177) and (sunpos.ra.degree < 183)
    if is_spring_equinox or is_fall_equinox:
        step_deg = 0.01
    else:
        step_deg = 0.1

    nsteps = int(np.round(360/step_deg))
    for i in range(nsteps):
        phi_now = phis[-1]
        theta_now = thetas[-1]
        
        # first of all, we assume that phi is uniformly scanned, but this is not the case
        # so we use this as an initial guess, and adjust thhe step, such that
        # the angle between consecutive points are the same (~1 degree, 360 points)
        
        step_rad = step_deg/180*np.pi
        phi_trial = phi_now + step_rad
        theta_trial = get_theta_given_phi(a, b, c, phi_trial)
        cossep_trial = cos_in_sphere(theta_now, phi_now, theta_trial, phi_trial)
        sep_deg_trial = np.arccos(cossep_trial)/np.pi*180
#       Yuhan originally tried to make the step sizes equal, but it's not
#       necessary--just oversample.
#        while abs(sep_deg_trial-step_deg)>0.001:
#            print(abs(sep_deg_trial-step_deg))
#            step_deg /= sep_deg_trial
#            step_rad = step_deg/180*np.pi
#            phi_trial = phi_now + step_rad
#            theta_trial = get_theta_given_phi(a, b, c, phi_trial)
#            cossep_trial = cos_in_sphere(theta_now, phi_now, theta_trial, phi_trial)
#            sep_deg_trial = np.arccos(cossep_trial)/np.pi*180
        phis.append(phi_trial)
        thetas.append(theta_trial)
        sep_degs.append(sep_deg_trial)
    phis = np.array(phis)
    thetas = np.array(thetas)
    alphas = phis
    decs = np.pi/2 - thetas
    return alphas/np.pi*180, decs/np.pi*180

def get_srg_fields(tnow, fields):
    """Return ZTF field IDs overlapping the SRG scanning footprint tonight.

    First attempts to load a pre-computed pointing plan from disk via
    `load_SRG_pointing`; falls back to the analytical `SRG_pointing` model
    on any exception. Matches SRG pointing centres to ZTF fields within
    5.25° (field centre-to-corner distance for the 7.5° × 7.3° ZTF FoV).

    Parameters
    ----------
    tnow : astropy.time.Time
        Any time within the night of interest.
    fields : Fields
        Pre-initialised ZTF field grid with observability already computed.

    Returns
    -------
    list of int
        ZTF field IDs whose centres lie within 5.25° of at least one SRG
        pointing centre tonight.
    """

    try:
        srg_ra, srg_dec = load_SRG_pointing(tnow)
    except Exception as e:
        logger.exception(e)
        srg_ra, srg_dec = SRG_pointing(tnow)

    # we will use a simple nearest neighbor to avoid lots of nested
    # for loops, as in the original code by Yuhan.
    # with the finer sampling that should be sufficient.
    
    srg_sc = SkyCoord(srg_ra, srg_dec, unit=u.degree)
    
    cuts = fields.select_fields(dec_range=[-31.,90.], 
                           grid_id=0, 
                           observable_hours_range=[1.0, 24.])


    ztf_sc = fields._field_coords(cuts=cuts)

    idx, d2d, _ = srg_sc.match_to_catalog_sky(ztf_sc)

    # max separation should be field center to field corner
    # field of view is 7.5 deg x 7.3 deg
    # let's call it 5.25 deg
    max_sep = 5.25 * u.degree
    sep_constraint = d2d < max_sep

    df_idx = cuts.loc[cuts].index

    return fields.fields.loc[df_idx[idx[sep_constraint]]].index.unique().tolist()
