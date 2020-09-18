"""Special-case field selection functions.

All functions should take 
    the current astropy.time.Time 
    an ObsLogger instance
    a dictionary of other observing programs field_ids and field_selection_functions
and return a list of field_ids."""


def srg_selection(time, obs_log, other_program_fields):
    return [666]
