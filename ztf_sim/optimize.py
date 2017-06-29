

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

def request_set_optimize(df_metric, df, requests_allowed):
    """Identify which request sets to observe tonight.

    Decision variable is yes/no per request_id"""

    request_sets = df_metric.index.values
    slots = df_metric.columns.values

    nreqs = df['total_requests_tonight']

    # calculate n_usable slots
    df_usable = df_metric > 0.05
    n_usable = df_usable.sum(axis=1)
    n_usable.name = 'n_usable' 

    # sum df_metric down to one column
    metric_sum = (df_metric * df_usable).sum(axis=1)
    metric_sum.name = 'metric_sum'

    # merge additional useful info
    dfr = df[['program_id','total_requests_tonight']].join(n_usable).join(metric_sum)

    dfr['occupancy'] = dfr['total_requests_tonight']/dfr['n_usable']

    # Create an empty model
    m = Model('requests')

    # decision variable: yes or no for each request set
    yr_dict = m.addVars(df_metric.index,name='Yr',vtype=GRB.BINARY)
    yr_series = pd.Series(yr_dict,name='Yr')
    dfr = dfr.join(yr_series)

    m.setObjective(np.sum(dfr['Yr'] * dfr['metric_sum'] * dfr['occupancy']), 
        GRB.MAXIMIZE)


    # slot occupancy constraint: nreqs obs divided over nusable slots
    constr_avg_slot_occupancy = m.addConstrs(
        ((np.sum(df_usable[t]*dfr['occupancy'] * dfr['Yr'])
        <= max_exps_per_slot) for t in slots), "constr_avg_slot_occupancy")

    # program balance
    constr_balance = m.addConstrs(
        ((np.sum(dfr.loc[dfr['program_id'] == p, 'Yr'] * 
                 dfr.loc[dfr['program_id'] == p, 'total_requests_tonight']  )
        <= requests_allowed[p]) for p in requests_allowed.keys()), 
        "constr_balance")

    # Quick and dirty is okay!
    m.Params.TimeLimit = 30.

    m.update()

    m.optimize()

    if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
        raise ValueError("Optimization failure")


    # now get the decision variables
    dfr['Yr_val'] = dfr['Yr'].apply(lambda x: x.getAttr('x'))
    dfr['Yr_val'] = dfr['Yr_val'].astype(bool)

    return dfr.loc[dfr['Yr_val'],'program_id'].index


def slot_optimize(df_metric, df, requests_allowed):
    """Determine which slots to place the requests in.

    Decision variable is yes/no per request_id, slot"""

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
        <= max_exps_per_slot) for t in slots), "constr_nperslot")

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

