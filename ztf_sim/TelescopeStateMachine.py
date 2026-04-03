"""State machine for simulating observations."""

from transitions import Machine
from astropy.time import Time
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import logging
from .utils import *
from .constants import BASE_DIR, P48_loc, FILTER_IDS
from .constants import READOUT_TIME, EXPOSURE_TIME, FILTER_CHANGE_TIME, slew_time

class TelescopeStateMachine(Machine):
    """State machine modelling the P48 telescope during simulated observations.

    Built on the ``transitions`` library. Tracks the physical state of the
    telescope and advances simulated time as it slews, changes filters, and
    exposes. Also handles night/weather checks via `can_observe`.

    States
    ------
    ready
        Telescope is idle and can accept the next command.
    slewing
        Telescope axes are moving to a new field.
    changing_filters
        The filter wheel is rotating to a new filter.
    exposing
        Detector is integrating.
    cant_observe
        Telescope cannot observe (daytime or weather).

    Attributes
    ----------
    current_time : astropy.time.Time
        Simulated UTC clock.
    current_ha : astropy.units.Quantity
        Current hour angle of the telescope pointing.
    current_dec : astropy.units.Quantity
        Current declination of the telescope pointing.
    current_domeaz : astropy.units.Quantity
        Current dome azimuth.
    current_filter_id : int
        Currently mounted filter (1 = g, 2 = r, 3 = i).
    current_zenith_seeing : astropy.units.Quantity
        Current zenith seeing FWHM.
    historical_observability_year : int or None
        PTF historical weather year used for weather simulation (2009–2015),
        or ``None`` for perfect-weather runs.
    """

    def __init__(self, current_time=Time('2018-01-01', scale='utc',
                                         location=P48_loc),
                 current_ha=0. * u.deg, current_dec=33.36 * u.deg,
                 current_domeaz=180. * u.deg,
                 current_filter_id=2, filters=FILTER_IDS,
                 current_zenith_seeing=2.0 * u.arcsec,
                 target_skycoord=None,
                 historical_observability_year=2015):
        """Initialise the telescope state machine.

        Parameters
        ----------
        current_time : astropy.time.Time, optional
            Starting simulation time. Default is 2018-01-01 UTC at P48.
        current_ha : astropy.units.Quantity, optional
            Starting hour angle. Default is 0 deg (on meridian).
        current_dec : astropy.units.Quantity, optional
            Starting declination. Default is 33.36 deg (Palomar latitude).
        current_domeaz : astropy.units.Quantity, optional
            Starting dome azimuth. Default is 180 deg (south).
        current_filter_id : int, optional
            Starting filter. Default is 2 (r-band).
        filters : list of int, optional
            All valid filter IDs. Default is ``FILTER_IDS``.
        current_zenith_seeing : astropy.units.Quantity, optional
            Starting zenith seeing FWHM. Default is 2.0 arcsec.
        target_skycoord : astropy.coordinates.SkyCoord or None, optional
            Sky coordinate of the current target. Default is ``None``.
        historical_observability_year : int or None, optional
            PTF observing year (2009–2015) to use for weather simulation. Set
            to ``None`` for perfect-weather (twilight-only) simulations.
            Default is 2015.
        """

        # Define some states.
        states = ['ready', 'cant_observe',
                  'slewing', 'changing_filters', 'exposing']

        # define the transitions

        transitions = [
            {'trigger': 'start_slew', 'source': 'ready', 'dest': 'slewing',
                'after': ['process_slew', 'stop_slew'],
                'conditions': 'slew_allowed'},
            {'trigger': 'stop_slew', 'source': 'slewing', 'dest': 'ready'},
            # for now do not require filter changes to include a slew....
            {'trigger': 'start_filter_change', 'source': 'ready',
                'dest': 'changing_filters',
                'after': ['process_filter_change', 'stop_filter_change']},
            {'trigger': 'stop_filter_change', 'source': 'changing_filters',
                'dest': 'ready'},
            {'trigger': 'start_exposing', 'source': 'ready', 'dest': 'exposing',
                'after': ['process_exposure', 'stop_exposing']},
            {'trigger': 'stop_exposing', 'source': 'exposing', 'dest': 'ready'},
            # I would like to automatically set the cant_observe state from
            # start_exposing, but that doesn't seem to work.
            {'trigger': 'check_if_ready', 'source': ['ready', 'cant_observe'],
                'dest': 'ready', 'conditions': 'can_observe'},
            {'trigger': 'set_cant_observe', 'source': '*',
                'dest': 'cant_observe'}
        ]

        # Initialize the state machine.  syntax from
        # https://github.com/tyarkoni/transitions
        Machine.__init__(self, states=states,
                         transitions=transitions,
                         initial='ready')

        self.current_time = current_time
        self.current_ha = current_ha
        self.current_dec = current_dec
        self.current_domeaz = current_domeaz
        self.current_filter_id = current_filter_id
        self.filters = filters
        self.current_zenith_seeing = current_zenith_seeing
        self.target_skycoord = target_skycoord

        # historical observability
        self.historical_observability_year = historical_observability_year
        self.observability = PTFObservabilityDB()

        self.logger = logging.getLogger(__name__)
        #self.logger = logging.getLogger('transitions')

    def current_state_dict(self):
        """Return a snapshot of the current telescope state as a dictionary.

        Returns
        -------
        dict
            Keys: ``'current_time'`` (astropy.time.Time),
            ``'current_ha'`` (astropy.units.Quantity),
            ``'current_dec'`` (astropy.units.Quantity),
            ``'current_domeaz'`` (astropy.units.Quantity),
            ``'current_filter_id'`` (int),
            ``'current_zenith_seeing'`` (astropy.units.Quantity),
            ``'filters'`` (list of int),
            ``'target_skycoord'`` (astropy.coordinates.SkyCoord or None).
        """
        return {'current_time': self.current_time,
                'current_ha': self.current_ha,
                'current_dec': self.current_dec,
                'current_domeaz': self.current_domeaz,
                'current_filter_id': self.current_filter_id,
                'current_zenith_seeing': self.current_zenith_seeing,
                'filters': self.filters,
                'target_skycoord': self.target_skycoord}

    def can_observe(self):
        """Check whether the telescope can currently observe.

        First tests the 12-degree evening/morning twilight constraint. If the
        Sun is above −12°, fast-forwards ``current_time`` to the next 12-degree
        evening twilight. If a ``historical_observability_year`` is set, also
        consults the PTF weather database; on a weathered-out block,
        fast-forwards to the end of the current block.

        Returns
        -------
        bool
            ``True`` if the telescope can observe at ``current_time``.
        """
        self.logger.info(self.current_time.iso)

        # start by checking for 12 degree twilight
        if coord.get_sun(self.current_time).transform_to(
                coord.AltAz(obstime=self.current_time,
                            location=P48_loc)).alt.is_within_bounds(
                upper=-12. * u.deg):
            if self.historical_observability_year is None:
                # don't use weather, just use 12 degree twilight
                return True
            else:
                is_observable = self.observability.check_historical_observability(
                    self.current_time, year=self.historical_observability_year)
                if not is_observable:
                    # optimization: fast-forward to start of next block
                    block_now = block_index(self.current_time)
                    block_end_time = block_index_to_time(block_now,
                        self.current_time, where='end')[0]
                    self.logger.info('Weathered out.  Fast forwarding to end of this block: {}'.format(
                        block_end_time.iso))
                    self.current_time = block_end_time

                return is_observable
        else:
            # daytime
            # optimization: fast-forward to sunset
            next_twilight = next_12deg_evening_twilight(self.current_time)
            self.logger.info('Fast forwarding to 12 deg twilight: {}'.format(
                next_twilight.iso))
            self.current_time = next_twilight
            return False

    def slew_allowed(self, target_skycoord):
        """Check whether a slew to *target_skycoord* is within telescope limits.

        Parameters
        ----------
        target_skycoord : astropy.coordinates.SkyCoord
            Target sky coordinate.

        Returns
        -------
        bool
            ``True`` if the target is reachable. ``False`` if the target
            altitude is below 10° or the declination is outside [−35°, +90°].
        """

        if (skycoord_to_altaz(target_skycoord, self.current_time).alt
            < (10. * u.deg)):
            return False

        if ((target_skycoord.dec < -35. * u.deg) or
                (target_skycoord.dec > 90. * u.deg)):
            return False
        return True

    def process_slew(self, target_skycoord,
                     readout_time=READOUT_TIME):
        """Advance time and update pointing after a slew.

        Evaluates the HA, Dec, and dome axes independently using the P48 slew
        model and takes the maximum. Advances ``current_time`` by
        ``max(axis_slew_times, readout_time)``. Updates ``current_ha``,
        ``current_dec``, and ``current_domeaz`` to the post-slew values (HA
        and dome Az are recomputed after the slew completes, accounting for
        sky rotation during the slew).

        Parameters
        ----------
        target_skycoord : astropy.coordinates.SkyCoord
            Sky position to slew to.
        readout_time : astropy.units.Quantity, optional
            Readout time of the previous exposure, which runs concurrently
            with the slew. Sets the minimum inter-exposure gap. Default is
            ``READOUT_TIME``.
        """
        # if readout_time is nonzero, assume we are reading during the slew,
        # which sets the lower limit for the time between exposures.

        self.target_skycoord = target_skycoord

        target_ha = RA_to_HA(self.target_skycoord.ra, self.current_time)
        target_domeaz = skycoord_to_altaz(self.target_skycoord,
                                          self.current_time).az
        target_dec = target_skycoord.dec

        # calculate time required to slew
        # duplicates codes in fields.py--consider refactoring
        axis_slew_times = [READOUT_TIME]
        for axis in ['ha', 'dec', 'domeaz']:
            dangle = np.abs(eval("target_{}".format(axis)) -
                            eval("self.current_{}".format(axis)))
            angle = np.where(dangle < (360. * u.deg - dangle), dangle,
                             360. * u.deg - dangle)
            axis_slew_times.append(slew_time(axis[:4], angle))

        net_slew_time = np.max([st.value for st in axis_slew_times]) *\
            axis_slew_times[0].unit

        # update the time
        self.current_time += net_slew_time
        # small deviation here: ha, az of target ra shifts (usually!)
        # modestly during slew,
        # so store the value after the slew is complete.

        target_ha = RA_to_HA(self.target_skycoord.ra, self.current_time)
        target_domeaz = skycoord_to_altaz(self.target_skycoord,
                                          self.current_time).az
        self.current_ha = target_ha
        self.current_dec = self.target_skycoord.dec
        self.current_domeaz = target_domeaz

    def process_filter_change(self, target_filter_id,
                              filter_change_time=FILTER_CHANGE_TIME):
        """Execute a filter change, advancing the clock if the filter differs.

        If ``current_filter_id`` already equals *target_filter_id*, no time
        is added.

        Parameters
        ----------
        target_filter_id : int
            Filter to switch to (1 = g, 2 = r, 3 = i).
        filter_change_time : astropy.units.Quantity, optional
            Time required to change the filter. Default is
            ``FILTER_CHANGE_TIME`` (135 s).
        """
        if self.current_filter_id != target_filter_id:
            self.current_filter_id = target_filter_id
            self.current_time += filter_change_time

    def process_exposure(self, exposure_time):
        """Advance time by *exposure_time* and update HA/dome for tracking.

        After the exposure completes, recalculates the current HA and dome
        azimuth of ``target_skycoord`` at the new ``current_time`` to account
        for sky rotation during the integration.

        Parameters
        ----------
        exposure_time : astropy.units.Quantity
            Duration of the exposure.
        """
        # annoyingly, transitions doesn't let me modify object
        # variables in the trigger functions themselves
        self.current_time += exposure_time
        # update ha and domeaz for tracking during the exposure
        target_ha = RA_to_HA(self.target_skycoord.ra, self.current_time)
        target_domeaz = skycoord_to_altaz(self.target_skycoord,
                                          self.current_time).az
        self.current_ha = target_ha
        self.current_domeaz = target_domeaz

    def wait(self, wait_time=EXPOSURE_TIME):
        """Advance the simulation clock without taking an observation.

        Used when the queue is empty, a slew fails, or the telescope is in
        ``cant_observe`` state.

        Parameters
        ----------
        wait_time : astropy.units.Quantity, optional
            Duration to wait. Default is ``EXPOSURE_TIME`` (30 s).
        """
        self.current_time += wait_time


class PTFObservabilityDB(object):
    """Historical PTF weather database for realistic weather simulation.

    Loads PTF observing records binned into 20-minute time blocks
    (``data/weather_blocks_20min.db``) and checks whether PTF was actively
    observing at the equivalent time in a specified historical year.

    Attributes
    ----------
    df : pandas.DataFrame
        Block-level PTF observation counts, indexed by ``(year, block)``.
    """

    def __init__(self):
        """Load the PTF weather database into memory.

        Reads ``data/weather_blocks_20min.db`` via `df_read_from_sqlite` and
        sets the index to ``(year, block)``.
        """
        df = df_read_from_sqlite('weather_blocks')
        self.df = df.set_index(['year', 'block'])

    def check_historical_observability(self, time, year=2015, nobs_min=5):
        """Given a (possibly future) UTC time, look up whether PTF 
        was observing at that time in specified year.

        Parameters
        ----------
        time : scalar astropy Time object
            UTC Time
        year : int [2009 -- 2015]
            year to check PTF historical observing
        nobs_min : int (default = 3)
            minimum number of observations per block to count as observable

        Returns
        -------
        bool
            ``True`` if the number of PTF observations in the matching block
            is greater than or equal to *nobs_min*; ``False`` if the block is
            absent (no observations) or below the threshold.
        """

        assert((year >= 2009) and (year <= 2015))
        assert((nobs_min > 0))

        block = block_index(time, time_block_size=TIME_BLOCK_SIZE)

        try:
            return self.df.loc[(year, block[0])].values[0] >= nobs_min
        except KeyError:
            # blocks are unfilled if there are no observations
            return False
