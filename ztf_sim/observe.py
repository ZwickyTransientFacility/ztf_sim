from ZTFStateMachine import ZTFStateMachine
from astropy.time import Time
import astropy.units as u


def observe():

    tel = ZTFStateMachine()

    while tel.current_time < Time('2018-01-02',scale='utc'):
        
        if tel.check_if_ready():
            # get coords
            # TODO: filter change
            if not tel.start_slew(0*u.deg,0*u.deg,0*u.deg):
                tel.set_cant_observe()
                tel.wait()
                continue
            if not tel.start_exposing():
                tel.set_cant_observe()
                tel.wait()
                continue
        else:
            tel.set_cant_observe()
            tel.wait()
            

