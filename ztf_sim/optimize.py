

from gurobipy import *
import numpy as np
import shelve
import astropy.units as u
import pandas as pd
from constants import TIME_BLOCK_SIZE, EXPOSURE_TIME, READOUT_TIME

#s = shelve.open('tmp_vars.shelf',flag='r')
#df_metric = s['block_slot_metric']
#df = s['df']

requests_allowed = {1: 548, 2: 548, 3: 274}

max_exps_per_slot = np.ceil((TIME_BLOCK_SIZE / 
                (EXPOSURE_TIME + READOUT_TIME)).to(
                u.dimensionless_unscaled).value).astype(int)

def slot_optimize(df_metric, df, requests_allowed):

    request_sets = df_metric.index.values
    slots = df_metric.columns.values

    nreqs = df['total_requests_tonight']

    # flatten the metric dataframe to make it easier to work with 
    df_metric['request_id'] = df_metric.index

    # make a "tidy" dataframe with one row per (request, slot [, filter])
    dft = pd.melt(df_metric,id_vars='request_id',var_name='slot',value_name='metric')

    dft = pd.merge(dft,df[['program_id','total_requests_tonight']],
        left_on='request_id',right_index=True)


    # Create an empty model
    m = Model('slots')

    yrt_dict = m.addVars(dft.index,name='Yrt',vtype=GRB.BINARY)
    yrt_series = pd.Series(yrt_dict,name='Yrt')
    dft = dft.join(yrt_series)

    def f_n(dft):
        grp = dft.groupby('request_id')
        f = (grp['Yrt'].agg(np.sum)/nreqs)
        f.name = 'f'
        return dft.join(f,on='request_id')['f']

    # this makes it quadratic, and slow
    #m.setObjective(np.sum(dft['Yrt'] * dft['metric'] * f_n(dft)), 
    #    GRB.MAXIMIZE)
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
        == max_exps_per_slot) for t in slots), "constr_nperslot")

    # program balance
    constr_balance = m.addConstrs(
        ((np.sum(dft.loc[dft['program_id'] == p, 'Yrt'])
        <= requests_allowed[p]) for p in requests_allowed.keys()), 
        "constr_balance")

    m.update()

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        raise ValueError("Optimization failure")


    # now get the decision variables
    dft['Yrt_val'] = dft['Yrt'].apply(lambda x: x.getAttr('x'))
    dft['Yrt_val'] = dft['Yrt_val'].astype(bool)

    return dft.loc[dft['Yrt_val'],['slot','request_id']]

