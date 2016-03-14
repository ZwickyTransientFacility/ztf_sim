
"""Functions defining cadence windows.  Return True if the field can
be observed at the supplied time."""



def no_cadence(*args):
    """No cadence requirement--can be observed at any time."""
    return True

def time_since_last_obs(time, pars, obshist):
    """

    window_start : time delta (astropy quantity)
        time since last observation after which observations can begin
    window_stop : time delta (astropy quantity)
        time since last observation after which observations stop.
        Set to np.inf * u.days to have an open-ended window.
    prev_filter : {'same','any', 'other', [specific filterkey]}
        which filter to use to determine time of last observation
    """
    
    pass
