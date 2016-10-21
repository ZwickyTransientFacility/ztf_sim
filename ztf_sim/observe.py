from ZTFStateMachine import ZTFStateMachine
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from QueueManager import GreedyQueueManager
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

run_name = 'tmp'


def observe(run_name=run_name, start_time='2016-03-20 02:30:00',
            weather_year=None, survey_duration=12. * u.hour,
            profile=profile):

    if profile:
        if survey_duration > 1. * u.day:
            print("Don't profile long runs: 25% overhead")
            profile = False
        else:
            profiler = Profiler()
            profiler.start()

    survey_start_time = Time(start_time, scale='utc', location=P48_loc)

    tel = ZTFStateMachine(
        current_time=survey_start_time,
        historical_observability_year=weather_year,
        logfile='../sims/log_{}'.format(run_name))

    # set up QueueManager
    Q = GreedyQueueManager()

    # set up Observing Programs
    CollabOP = CollaborationObservingProgram(
        Q.fields.select_field_ids(dec_range=[-30, 90], abs_b_range=[20, 90],
                                  grid_id=0))
    Q.add_observing_program(CollabOP)
    MSIPOP = MSIPObservingProgram(
        Q.fields.select_field_ids(dec_range=[-30, 90], grid_id=0))
    #MSIPOP.observing_time_fraction = 1.0
    Q.add_observing_program(MSIPOP)
    CaltechOP = CaltechObservingProgram(
        Q.fields.select_field_ids(dec_range=[-30, 90], grid_id=0))
    Q.add_observing_program(CaltechOP)
    # initialize nightly field requests (Tom Barlow function)
    Q.assign_nightly_requests(tel.current_state_dict())

    # initialize sqlite history
    log = ObsLogger(run_name, tel.current_time)
    log.create_pointing_log(clobber=True)

    current_night_mjd = np.floor(tel.current_time.mjd)

    while tel.current_time < (survey_start_time + survey_duration):

        # check if it is a new night and reload queue with new requests
        if np.floor(tel.current_time.mjd) > current_night_mjd:
            Q.assign_nightly_requests(tel.current_state_dict())
            current_night_mjd = np.floor(tel.current_time.mjd)

        if tel.check_if_ready():
            current_state = tel.current_state_dict()
            # get coords
            next_obs = Q.next_obs(current_state)

            # try to change filters, if needed
            if next_obs['target_filter_id'] != current_state['current_filter_id']:
                if not tel.start_filter_change(next_obs['target_filter_id']):
                    # TODO: log the failure
                    log.prev_obs = None
                    tel.wait()
                    continue

            # try to slew to the next target
            if not tel.start_slew(coord.SkyCoord(next_obs['target_ra'] * u.deg,
                                                 next_obs['target_dec'] * u.deg)):
                tel.set_cant_observe()
                # TODO: log the failure
                # "missed history": http://ops2.lsst.org/docs/current/architecture.html#output-tables
                log.prev_obs = None
                tel.wait()
                continue

            # try to expose
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
                # c) remove completed request_id from the pool and the queue
                Q.remove_requests(next_obs['request_id'])
        else:
            tel.set_cant_observe()
            tel.wait()

    if profile:
        profiler.stop()
        print profiler.output_text(unicode=True, color=True)
        with open('../sims/profile_{}'.format(run_name), 'w') as f:
            f.write(profiler.output_text())
