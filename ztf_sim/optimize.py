

from gurobipy import *
import numpy as np
import shelve
from constants import TIME_BLOCK_SIZE, EXPOSURE_TIME, READOUT_TIME

s = shelve.open('tmp_vars.shelf',flag='r')
df_metric = s['block_slot_metric']
df = s['df']
nreqs = df['total_requests_tonight']

max_exps_per_slot = np.ceil((TIME_BLOCK_SIZE / 
                (EXPOSURE_TIME + READOUT_TIME)).to(
                u.dimensionless_unscaled).value).astype(int)

def slot_optimize():

    # Create an empty model
    m = Model('slots')

    request_sets = df_metric.index.values
    slots = df_metric.columns.values
    # add a variable f
    
    Yrt = m.addVars(request_sets, slots, vtype=GRB.BINARY, obj=df_metric.values,
        name="Yrt")
    # adding the objective function here is useful iff we're just doing a
    # direct sum, with no f(n)

    m.setObjective(Yrt.sum(), GRB.MAXIMIZE)

    # no more than nreqs slots  assigned per request set
    constr_nreqs = m.addConstrs(
        (Yrt.sum(r,'*') <= nreqs[r] for r in request_sets), "constr_nreqs")


    # no more than nobs per slot
    # TODO: tune this so it's closer to the average performance so we don't
    # drop too many...
    constr_nperslot = m.addConstrs(
        (Yrt.sum('*',t) <= max_exps_per_slot for t in slots), "constr_nperslot")

    m.update()

    m.optimize()
