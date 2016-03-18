from ZTFStateMachine import ZTFStateMachine
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from queue import QueueManager


# TODO: set up configuration system so we can easily run (and distinguish)
# sims with various inputs.  tag with commit hash!

def observe():

    tel = ZTFStateMachine()

    # set up QueueManager with field requests (Tom Barlow function)
    # reload each night?
    Q = GreedyQueueManager()

    # TODO: initialize sqlite history

    while tel.current_time < Time('2018-01-02',scale='utc'):
        
        if tel.check_if_ready():
            current_state = tel.current_state()
            # get coords
            next_obs = Q.next_obs(current_state)

            # TODO: filter change, if needed
            
            if not tel.start_slew(
                    coord.Skycoord(next_obs['ra'],next_obs['dec'])):
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
                Q.remove_requests(next_obs['request_id'])
        else:
            tel.set_cant_observe()
            tel.wait()
            

