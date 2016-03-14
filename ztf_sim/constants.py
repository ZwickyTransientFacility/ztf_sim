
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
import astropy.units as u

P48_loc = coords.EarthLocation(lat=coords.Latitude('33d21m26.35s'),
    lon=coords.Longitude('-116d51m32.04s'), height=1707.)

P48_slew_pars = {
    'ha':{'coord':'ra','accel':0.27*u.deg*u.second**(-2.),
        'decel':0.27*u.deg*u.second**(-2.),'vmax':1.6*u.deg/u.second},
    'dec':{'coord':'dec','accel':0.2*u.deg*u.second**(-2.),
        'decel':0.15*u.deg*u.second**(-2.),'vmax':1.2*u.deg/u.second},
    'dome':{'coord':'az','accel':0.17*u.deg*u.second**(-2.),
        'decel':0.6*u.deg*u.second**(-2.),'vmax':3.*u.deg/u.second}}

EXPOSURE_TIME = 30.*u.second
READOUT_TIME = 10.*u.second
FILTER_CHANGE_TIME = 90.*u.second
SETTLE_TIME = 2.*u.second

TIME_BLOCK_SIZE = 20.*u.min

program_name_to_id = {'collaboration':1, 'MSIP':2, 'Caltech':3}
program_names = program_name_to_id.keys()
program_id_to_name = {v: k for k, v in program_name_to_id.items()}
program_ids = program_id_to_name.keys()

filter_name_to_id = {'g':1, 'r':2}
filter_names = filter_name_to_id.keys()
filter_id_to_name = {v: k for k, v in filter_name_to_id.items()}
filter_ids = filter_id_to_name.keys()


def slew_time(axis, angle):
    vmax = P48_slew_pars[axis]['vmax']
    acc = P48_slew_pars[axis]['accel']
    dec = P48_slew_pars[axis]['decel']

    t_acc=vmax/acc
    t_dec=vmax/dec
    # TODO: if needed, include conditional for scalars
#    if 0.5*vmax*(t_acc+t_dec)>=angle:
#        return np.sqrt(2*angle*(1./acc+1./dec))
#    else:
#        return 0.5*(2.*angle/vmax+t_acc+t_dec)
    slew_time = 0.5*(2.*angle/vmax+t_acc+t_dec)
    w = 0.5*vmax*(t_acc+t_dec) >= angle
    slew_time[w] = np.sqrt(2*angle[w]*(1./acc+1./dec))
    return slew_time + SETTLE_TIME
