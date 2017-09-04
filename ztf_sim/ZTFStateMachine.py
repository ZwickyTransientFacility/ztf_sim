from __future__ import absolute_import
from builtins import object
from transitions import Machine
from astropy.time import Time
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from utils import *
from constants import *
import logging
from transitions import logger


class ZTFStateMachine(Machine):

    def __init__(self, current_time=Time('2018-01-01', scale='utc',
                                         location=P48_loc),
                 current_ha=0. * u.deg, current_dec=33.36 * u.deg,
                 current_domeaz=180. * u.deg,
                 current_filter_id=2, filters=FILTER_IDS,
                 current_zenith_seeing=2.0 * u.arcsec,
                 target_skycoord=None,
                 logfile='../sims/log_ztf_sim',
                 historical_observability_year=2015):

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

        # logging.  wipe out existing log.
        fh = logging.FileHandler(logfile, mode='w')
        fh.setLevel(logging.INFO)
        self.logger = logger
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(fh)

    def current_state_dict(self):
        """Return current state parameters in a dictionary"""
        return {'current_time': self.current_time,
                'current_ha': self.current_ha,
                'current_dec': self.current_dec,
                'current_domeaz': self.current_domeaz,
                'current_filter_id': self.current_filter_id,
                'current_zenith_seeing': self.current_zenith_seeing,
                'filters': self.filters,
                'target_skycoord': self.target_skycoord}

    def can_observe(self):
        """Check for night and weather"""
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
        """Check that slew is within allowed limits"""

        if (skycoord_to_altaz(target_skycoord, self.current_time).alt
            < (10. * u.deg)):
            return False

        if ((target_skycoord.dec < -35. * u.deg) or
                (target_skycoord.dec > 90. * u.deg)):
            return False
        return True

    def process_slew(self, target_skycoord,
                     readout_time=READOUT_TIME):
        # if readout_time is nonzero, assume we are reading during the slew,
        # which sets the lower limit for the time between exposures.

        self.target_skycoord = target_skycoord

        target_ha = RA_to_HA(self.target_skycoord.ra, self.current_time)
        target_domeaz = skycoord_to_altaz(self.target_skycoord,
                                          self.current_time).az
        target_dec = target_skycoord.dec

        # calculate time required to slew
        # TODO: duplicates codes in fields.py--consider refactoring
        axis_slew_times = [READOUT_TIME]
        for axis in ['ha', 'dec', 'domeaz']:
            dangle = np.abs(eval("target_{}".format(axis)) -
                            eval("self.current_{}".format(axis)))
            angle = np.where(dangle < (360. * u.deg - dangle), dangle,
                             360. * u.deg - dangle)
            axis_slew_times.append(slew_time(axis[:4], angle * u.deg))

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
        if self.current_filter_id != target_filter_id:
            self.current_filter_id = target_filter_id
            self.current_time += filter_change_time
        # TODO: put in actual treatment of filter change (e.g., slew to stow
        # position)

    def process_exposure(self, exposure_time=EXPOSURE_TIME):
        # annoyingly, transitions doesn't let me modify object
        # variables in the trigger functions themselves
        self.current_time += exposure_time
        # update ha and domeaz for tracking during the exposure
        target_ha = RA_to_HA(self.target_skycoord.ra, self.current_time)
        target_domeaz = skycoord_to_altaz(self.target_skycoord,
                                          self.current_time).az
        self.current_ha = target_ha
        self.current_domeaz = target_domeaz

        # TODO: put in logic that near-zenith pointings could create
        # long and slow dome slews during the exposure

    def wait(self, wait_time=EXPOSURE_TIME):
        self.current_time += wait_time


class PTFObservabilityDB(object):

    def __init__(self):
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
        ----------
        boolean if nobs > nobs_min
        """

        assert((year >= 2009) and (year <= 2015))
        assert((nobs_min > 0))

        block = block_index(time, time_block_size=TIME_BLOCK_SIZE)

        try:
            return self.df.loc[(year, block[0])].values[0] >= nobs_min
        except TypeError as KeyError:
            # blocks are unfilled if there are no observations
            return False
