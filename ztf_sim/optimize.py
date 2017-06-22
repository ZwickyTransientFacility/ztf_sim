

from gurobipy import *
import numpy as np
import shelve
import astropy.units as u
import pandas as pd
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
    # probably need to use tupledict's select and prod calls otherwise

    m.setObjective(Yrt.sum(), GRB.MAXIMIZE)

    # no more than nreqs slots  assigned per request set
    constr_nreqs = m.addConstrs(
        (Yrt.sum(r,'*') <= nreqs[r] for r in request_sets), "constr_nreqs")

    # no more than nobs per slot
    # TODO: tune this so it's closer to the average performance so we don't
    # drop too many...
    constr_nperslot = m.addConstrs(
        (Yrt.sum('*',t) <= max_exps_per_slot for t in slots), "constr_nperslot")

 .update()

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        raise ValueError("Optimization failure")

    # an ugly constructor
    df_Yrt = (df_metric * 0.).astype(bool)

    # is there a way to get the decision variable out of their tupledict
    # without a double loop????
    # using loc to output this is super slow, need another approach
    for k, v in Yrt.iteritems():
        df_Yrt.loc[k] = v.X

    df_Yrt = df_Yrt.astype(bool)

def pandas_optimize():

    request_sets = df_metric.index.values
    slots = df_metric.columns.values

    # flatten the metric dataframe to make it easier to work with 
    df_metric['request_id'] = df_metric.index

    # make a "tidy" dataframe with one row per (request, slot [, filter])
    dft = pd.melt(df_metric,id_vars='request_id',var_name='slot',value_name='metric')

    # Create an empty model
    m = Model('slots')

    # add a variable f
    
    # this is slow because of the iterrorws
    #yrt_list = [m.addVar(name="Yrt[{},{}]".format(row.request_id,row.slot),
    #            vtype=GRB.BINARY) for idxi, row in dft.iterrows()]
    #dft['Yrt'] = yrt_list

    # this doesn't work right
    #dft.assign(var=lambda x: m.addVar(name="Yrt[{},{}]".format(x.request_id,x.slot)))

    yrt_dict = m.addVars(dft.index,name='Yrt',vtype=GRB.BINARY)
    yrt_series = pd.Series(yrt_dict,name='Yrt')
    dft = dft.join(yrt_series)

    m.setObjective(np.sum(dft['Yrt'] * dft['metric']), 
        GRB.MAXIMIZE)

    # no more than nreqs slots  assigned per request set
    constr_nreqs = m.addConstrs(
        ((np.sum(dft.loc[dft['request_id'] == r, 'Yrt'])
        <= nreqs[r]) for r in request_sets), "constr_nreqs")


    # no more than nobs per slot
    # TODO: tune this so it's closer to the average performance so we don't
    # drop too many...
    constr_nperslot = m.addConstrs(
        ((np.sum(dft.loc[dft['slot'] == t, 'Yrt'])
        <= max_exps_per_slot) for t in slots), "constr_nperslot")

    m.update()

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        raise ValueError("Optimization failure")


    # now get the decision variables
    dft['Yrt_val'] = dft['Yrt'].apply(lambda x: x.getAttr('x'))
    dft['Yrt_val'] = dft['Yrt_val'].astype(bool)
