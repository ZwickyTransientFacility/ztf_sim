
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
import astropy.units as u
import astroplan

P48_loc = coords.EarthLocation(lat=coords.Latitude('33d21m26.35s'),
                               lon=coords.Longitude('-116d51m32.04s'),
                               height=1707.)

# use UTC only
P48_Observer = astroplan.Observer(location=P48_loc)

# HA and Dec from http://www.oir.caltech.edu/twiki_oir/bin/view/Palomar/ZTF/TelescopeSpecifications v5
# Dome estimate from Jeff Z email, 9/21/15
# Requirements/goals info from Jeff Z email, 12/12/16
P48_slew_pars = {
    'ha': {'coord': 'ra', 'accel': 0.27 * u.deg * u.second**(-2.),
           'decel': 0.27 * u.deg * u.second**(-2.),
           'vmax': 1.18 * u.deg / u.second},
    'dec': {'coord': 'dec', 'accel': 0.41 * u.deg * u.second**(-2.),
            'decel': 0.41 * u.deg * u.second**(-2.),
            'vmax': 1.50 * u.deg / u.second},
    'dome': {'coord': 'az', 'accel': 0.5 * u.deg * u.second**(-2.),
             'decel': 0.5 * u.deg * u.second**(-2.),
             'vmax': 5. * u.deg / u.second}}
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

EXPOSURE_TIME = 30. * u.second
READOUT_TIME = 10. * u.second
FILTER_CHANGE_TIME = 90. * u.second
SETTLE_TIME = 2. * u.second

MAX_AIRMASS = 2.5

TIME_BLOCK_SIZE = 20. * u.min

PROGRAM_NAME_TO_ID = {'collaboration': 1, 'MSIP': 2, 'Caltech': 3}
PROGRAM_NAMES = PROGRAM_NAME_TO_ID.keys()
PROGRAM_ID_TO_NAME = {v: k for k, v in PROGRAM_NAME_TO_ID.items()}
PROGRAM_IDS = PROGRAM_ID_TO_NAME.keys()

PROGRAM_BLOCK_SEQUENCE = [1, 2, 1, 2, 3]
LEN_BLOCK_SEQUENCE = len(PROGRAM_BLOCK_SEQUENCE)

FILTER_NAME_TO_ID = {'g': 1, 'r': 2, 'i': 3}
FILTER_NAMES = FILTER_NAME_TO_ID.keys()
FILTER_ID_TO_NAME = {v: k for k, v in FILTER_NAME_TO_ID.items()}
FILTER_IDS = FILTER_ID_TO_NAME.keys()

PIXEL_SCALE = 1.006  # arcsec/pixel


def slew_time(axis, angle):
    vmax = P48_slew_pars[axis]['vmax']
    acc = P48_slew_pars[axis]['accel']
    dec = P48_slew_pars[axis]['decel']

    t_acc = vmax / acc
    t_dec = vmax / dec
    # TODO: if needed, include conditional for scalars
#    if 0.5*vmax*(t_acc+t_dec)>=angle:
#        return np.sqrt(2*angle*(1./acc+1./dec))
#    else:
#        return 0.5*(2.*angle/vmax+t_acc+t_dec)
    slew_time = 0.5 * (2. * angle / vmax + t_acc + t_dec)
    w = 0.5 * vmax * (t_acc + t_dec) >= angle
    slew_time[w] = np.sqrt(2 * angle[w] * (1. / acc + 1. / dec))

    wnonzero = slew_time > 0
    slew_time[wnonzero] += SETTLE_TIME
    return slew_time 
