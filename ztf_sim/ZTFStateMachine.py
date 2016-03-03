from transitions import Machine
from astropy.time import Time
import numpy as np
import astropy.units as u
from utils import slew_time, READOUT_TIME, EXPOSURE_TIME

class ZTFStateMachine(object):

    # Define some states. 
    states = ['ready', 'cant_observe', 
            'slewing', 'changing_filters', 'exposing']

    # define the transitions
    
    transitions = [
        {'trigger': 'start_slew', 'source': 'ready', 'dest': 'slewing',
            'after': 'stop_slew', 'conditions': 'slew_allowed'},
        {'trigger': 'stop_slew', 'source': 'slewing', 'dest': 'ready'},
        # for now do not require filter changes to include a slew....
        {'trigger': 'start_filter_change', 'source': 'ready', 
            'dest': 'changing_filters', 'after': 'stop_filter_change'},
        {'trigger': 'stop_filter_change', 'source': 'changing_filters', 
            'dest': 'ready'},
        {'trigger': 'start_expose', 'source': 'ready', 'dest': 'exposing',
            'after': 'stop_expose', 'conditions': 'check_can_observe'},
        {'trigger': 'stop_expose', 'source': 'exposing', 'dest': 'ready'},
        {'trigger': 'wait', 'source': ['ready','cant_observe'], 
            'dest': 'ready'},
        {'trigger': 'set_cant_observe', 'source':'*', 'dest': 'cant_observe'}
        ]


    def __init__(self, current_time = Time('2018-01-01',scale='utc'), 
            current_ha = 0.*u.deg, current_dec = 0.*u.deg,
            current_domeaz = 0.*u.deg,
            current_filter = 'r', filters = ['r','g'],
            current_fieldid = None):

        self.current_time = current_time
        self.current_ha = current_ha
        self.current_dec = current_dec
        self.current_domeaz = current_domeaz
        self.current_filter = current_filter
        self.filters = filters
        self.current_fieldid = None

        # Initialize the state machine.  syntax from
        # https://github.com/tyarkoni/transitions
        self.machine = Machine(model=self, states=ZTFStateMachine.states, 
                initial='ready', transitions=ZTFStateMachine.transitions)

    def check_can_observe(self):
        """Check for night and weather"""
        # TODO: figure out if this is the right way to set this up...
        if False:
            self.set_cant_observe()
        return True

    def slew_allowed(self, target_ha, target_dec, target_domeaz):
        """Check that slew is within allowed limits"""
        if np.abs(target_ha) > 90.*u.deg:
            return False
        if target_dec < -40.*u.deg: # TODO: fix this value
            return False
        return True

    def start_slew(self, target_ha, target_dec, target_domeaz,
            readout_time = READOUT_TIME):
        # small deviation here: ha of target ra shifts modestly during slew
        # if readout_time is nonzero, assume we are reading during the slew,
        # which sets the lower limit for the time between exposures.
        
        # calculate time required to slew
        axis_slew_times = [READOUT_TIME]
	for axis in ['ha', 'dec', 'domeaz']:
            dangle = np.abs(eval("target_{}".format(axis)) -
                    eval("self.current_{}".format(axis)))
            angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
            axis_slew_times.append(slew_time(axis[:4], angle*u.deg))

        slew_time = np.max(axis_slew_times)

        # update the time
        self.current_time += slew_time
        self.current_ha = target_ha
        self.current_dec = target_dec
        self.current_domeaz = target_domeaz


    def start_exposing(self, exposure_time = EXPOSURE_TIME):
        self.current_time += exposure_time
        # TODO: update ra, dec, domeaz for shifts during exposure
        # TODO: put the logging here? (or as a separate callback)
        
    def wait(self, wait_time = EXPOSURE_TIME):
        self.current_time += wait_time
