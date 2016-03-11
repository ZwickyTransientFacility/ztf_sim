from ZTFStateMachine import ZTFStateMachine
from astropy.time import Time
import astropy.units as u
from queue import QueueManager


def observe():

    tel = ZTFStateMachine()

    # set up QueueManager with field requests
    Q = QueueManager()

    while tel.current_time < Time('2018-01-02',scale='utc'):
        
        if tel.check_if_ready():
            current_state = tel.current_state()
            # get coords
            next_obs = Q.next_obs(current_state)

            # TODO: filter change, if needed
            
            # where do I do the ra, dec to ha, dec, domeaz conversion?
            # I think it needs to be in the telescope state machine--need a
            # "tracking" state (or simulation).
            if not tel.start_slew(0*u.deg,0*u.deg,0*u.deg):
                tel.set_cant_observe()
                tel.wait()
                continue
            if not tel.start_exposing():
                tel.set_cant_observe()
                tel.wait()
                continue
            else:
                # exposure completed successfully.  now 
                # a) store exposure information in pointing history sqlite db
                # b) remove completed request_id
        else:
            tel.set_cant_observe()
            tel.wait()
            

