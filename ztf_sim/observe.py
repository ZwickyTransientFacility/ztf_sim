from ZTFStateMachine import ZTFStateMachine
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from queue import GreedyQueueManager
from ObsLogger import ObsLogger
from ObservingProgram import *
from constants import *

profile = False

if profile:
    try:
        from pyinstrument import Profiler
    except ImportError:
        print 'Error importing pyinstrument'
        profile = False

# TODO: set up configuration system so we can easily run (and distinguish)
# sims with various inputs.  tag with commit hash!
# or sub-tables of the db output...

run_name = 'test_run'


def observe(run_name=run_name):

    if profile:
        profiler = Profiler()
        profiler.start()

    tel = ZTFStateMachine(
        current_time=Time('2018-01-01 04:00:00', scale='utc',
                          location=P48_loc),
        # no weather
        historical_observability_year=None)

    # set up QueueManager
    Q = GreedyQueueManager()

    # set up Observing Programs
    #CollabOP = CollaborationObservingProgram()
    # Q.add_observing_program(CollabOP)
    MSIPOP = MSIPObservingProgram(
        Q.fields.select_field_ids(dec_range=[-30, 90], grid_id=0))
    MSIPOP.observing_time_fraction = 1.0
    Q.add_observing_program(MSIPOP)
    #CaltechOP = CaltechObservingProgram()
    # Q.add_observing_program(CaltechOP)

    # initialize nightly field requests (Tom Barlow function)
    Q.assign_nightly_requests(tel.current_state_dict())

    # temporary loading to test things
    # Q.rp.add_requests(1,
    #        Q.fields.fields[
    #            Q.fields.select_fields(dec_range=[-30,90])].index, 2,
    #        'no_cadence',{})

    # initialize sqlite history
    log = ObsLogger(run_name, tel.current_time)
    log.create_pointing_log(clobber=True)

    # while tel.current_time < Time('2018-01-02',scale='utc'):
    while tel.current_time < Time('2018-01-01 05:00:00', scale='utc'):

        # TODO: reload queue with new requests on update interval (nightly

        if tel.check_if_ready():
            current_state = tel.current_state_dict()
            # get coords
            next_obs = Q.next_obs(current_state)

            # TODO: filter change, if needed

            if not tel.start_slew(coord.SkyCoord(next_obs['target_ra'] * u.deg,
                                                 next_obs['target_dec'] * u.deg)):
                tel.set_cant_observe()
                # TODO: log the failure
                # "missed history": http://ops2.lsst.org/docs/current/architecture.html#output-tables
                log.prev_obs = None
                tel.wait()
                continue
            if not tel.start_exposing():
                tel.set_cant_observe()
                # TODO: log the failure
                log.prev_obs = None
                tel.wait()
                continue
            else:
                # exposure completed successfully.  now
                # a) store exposure information in pointing history sqlite db
                current_state = tel.current_state_dict()
                log.log_pointing(current_state, next_obs)
                # b) update Fields
                Q.fields.mark_field_observed(next_obs, current_state)
                # c) remove completed request_id
                Q.rp.remove_requests(next_obs['request_id'])
        else:
            tel.set_cant_observe()
            tel.wait()

    if profile:
        profiler.stop()
        print profiler.output_text(unicode=True, color=True)
