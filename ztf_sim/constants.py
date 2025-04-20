"""Constants."""

import os
import inspect
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
import astropy.units as u
import astroplan

BASE_DIR = os.path.dirname(os.path.abspath(inspect.getfile(
                inspect.currentframe()))) + '/'

# Palomar
# Defined by Eric Bellm
P48_loc = coords.EarthLocation(lat=coords.Latitude('33d21m26.2s'),
                               lon=coords.Longitude('-116d51m35.5s'),
                               height=1707.)
# From Steve, 2025-01-20
# La Silla
# Taken from https://en.wikipedia.org/wiki/La_Silla_Observatory
P48_loc = coords.EarthLocation(lat=coords.Latitude('-29d15m27s'),
                              lon=coords.Longitude('-70d44m15s'),
                              height=2400.)

# use UTC only
P48_Observer = astroplan.Observer(location=P48_loc)

# HA and Dec from http://www.oir.caltech.edu/twiki_oir/bin/view/Palomar/ZTF/TelescopeSpecifications v5
# Dome estimate from Jeff Z email, 9/21/15
# goals info from Jeff Z email, 12/12/16
# Ha/Dec from Telescope Drive Performance Assessment v1.2; dome estimate from
# Jeff Z. email, 9/27/17
P48_slew_pars = {
    'ha': {'coord': 'ra', 'accel': 0.4 * u.deg * u.second**(-2.),
           'decel': 0.4 * u.deg * u.second**(-2.),
           'vmax': 2.5 * u.deg / u.second},
    'dec': {'coord': 'dec', 'accel': 0.5 * u.deg * u.second**(-2.),
            'decel': 0.5 * u.deg * u.second**(-2.),
            'vmax': 3.0 * u.deg / u.second},
    'dome': {'coord': 'az', 'accel': 0.5 * u.deg * u.second**(-2.),
             'decel': 0.5 * u.deg * u.second**(-2.),
             'vmax': 3. * u.deg / u.second}}

P48_slew_pars_goal = {
    'ha': {'coord': 'ra', 'accel': 0.50 * u.deg * u.second**(-2.),
           'decel': 0.50 * u.deg * u.second**(-2.),
           'vmax': 3.00 * u.deg / u.second},
    'dec': {'coord': 'dec', 'accel': 0.5 * u.deg * u.second**(-2.),
            'decel': 0.5 * u.deg * u.second**(-2.),
            'vmax': 3.00 * u.deg / u.second},
    'dome': {'coord': 'az', 'accel': 0.5 * u.deg * u.second**(-2.),
             'decel': 0.5 * u.deg * u.second**(-2.),
             'vmax': 5. * u.deg / u.second}}

# From Steve, 2025-01-20
# LS4 Camera/telescope properties

# PTF configuration
# P48_slew_pars = {
#    'ha': {'coord': 'ra', 'accel': 0.27 * u.deg * u.second**(-2.),
#           'decel': 0.27 * u.deg * u.second**(-2.),
#           'vmax': 1.6 * u.deg / u.second},
#    'dec': {'coord': 'dec', 'accel': 0.2 * u.deg * u.second**(-2.),
#            'decel': 0.15 * u.deg * u.second**(-2.),
#            'vmax': 1.2 * u.deg / u.second},
#    'dome': {'coord': 'az', 'accel': 0.17 * u.deg * u.second**(-2.),
#             'decel': 0.6 * u.deg * u.second**(-2.),
#             'vmax': 3. * u.deg / u.second}}

# ZTF
EXPOSURE_TIME = 30. * u.second
READOUT_TIME = 8. * u.second
FILTER_CHANGE_TIME = 135. * u.second
SETTLE_TIME = 1. * u.second
MAX_AIRMASS = 2.5


# From Steve, 2025-01-20
# LS4
EXPOSURE_TIME = 45. * u.second                  # from 2023 CfP
READOUT_TIME  = 15. * u.second                  # Assuming read-out in two directions.
FILTER_CHANGE_TIME = 0. * u.second
SETTLE_TIME = 1. * u.second
MAX_AIRMASS = 2.0

TIME_BLOCK_SIZE = 30. * u.min

PROGRAM_NAME_TO_ID = {'engineering': 0, 
                      'MSIP':1, 'collaboration': 2, 'Caltech': 3}

#PROGRAM_NAME_TO_ID = {'engineering': 0, 
#                      'public':1, 'collaboration': 2, 'ToO': 3}

PROGRAM_NAMES = list(PROGRAM_NAME_TO_ID.keys())
PROGRAM_ID_TO_NAME = {v: k for k, v in list(PROGRAM_NAME_TO_ID.items())}
PROGRAM_IDS = list(PROGRAM_ID_TO_NAME.keys())

PROGRAM_BLOCK_SEQUENCE = [1, 2, 1, 2, 3]
LEN_BLOCK_SEQUENCE = len(PROGRAM_BLOCK_SEQUENCE)

FILTER_NAME_TO_ID = {'g': 1, 'r': 2, 'i': 3}
FILTER_NAMES = list(FILTER_NAME_TO_ID.keys())
FILTER_ID_TO_NAME = {v: k for k, v in list(FILTER_NAME_TO_ID.items())}
FILTER_IDS = list(FILTER_ID_TO_NAME.keys())

PIXEL_SCALE = 1.006  # arcsec/pixel


def slew_time(axis, angle):
    vmax = P48_slew_pars[axis]['vmax']
    acc = P48_slew_pars[axis]['accel']
    dec = P48_slew_pars[axis]['decel']

    t_acc = vmax / acc
    t_dec = vmax / dec
    slew_time = 0.5 * (2. * angle / vmax + t_acc + t_dec)
    w = 0.5 * vmax * (t_acc + t_dec) >= angle
    slew_time[w] = np.sqrt(2 * angle[w] * (1. / acc + 1. / dec))

    wnonzero = slew_time > 0
    slew_time[wnonzero] += SETTLE_TIME
    return slew_time 
