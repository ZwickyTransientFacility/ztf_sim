from __future__ import print_function
from __future__ import absolute_import

import os
from builtins import str
from builtins import range
from gurobipy import *
import numpy as np
import shelve
import astropy.units as u
import pandas as pd
from collections import defaultdict
from .constants import TIME_BLOCK_SIZE, EXPOSURE_TIME, READOUT_TIME

#s = shelve.open('tmp_vars.shelf',flag='r')
#df_metric = s['block_slot_metric']
#df = s['df']
#
#requests_allowed = {1: 548, 2: 548, 3: 274}

max_exps_per_slot = np.ceil((TIME_BLOCK_SIZE / 
                (EXPOSURE_TIME + READOUT_TIME)).to(
                u.dimensionless_unscaled).value).astype(int)

def request_set_optimize(df_metric, df, requests_allowed):
    """Identify which request sets to observe tonight.

    Decision variable is yes/no per request_id"""

    request_sets = df_metric.index.values
    slots = np.unique(df_metric.columns.get_level_values(0).values)
    filter_ids = np.unique(df_metric.columns.get_level_values(1).values)

    # can try working with it in 2D/3D, but may be easier just tidy up
    #idx = pd.IndexSlice
    #df_metric.loc[idx[:],idx[:,2]]
    # df_metric.unstack()
    
    # make a copy so I don't have downstream problems
    df_metric_local = df_metric.copy()
    
    df_metric_local['request_id'] = df_metric_local.index
    dft = pd.melt(df_metric_local,id_vars='request_id',
        var_name=['slot','metric_filter_id'],
        value_name='metric')
    # get n_reqs by fid
    n_reqs_cols = ['n_reqs_{}'.format(fid) for fid in filter_ids]
    n_reqs_cols.append('total_requests_tonight')
    dft = pd.merge(dft,df[n_reqs_cols],left_on='request_id',right_index=True)

    # nreqs_{fid} weighted sum over the filters
    for fid in filter_ids:
        wfid = dft['metric_filter_id'] == fid
        n_req_col = 'n_reqs_{}'.format(fid)
        dft.loc[wfid, 'metric'] *= (dft.loc[wfid, n_req_col] / 
            dft.loc[wfid, 'total_requests_tonight'])
    grprs = dft.groupby(['request_id','slot'])
    dfrs = grprs['metric'].agg(np.sum)

    # calculate n_usable slots
    grpr = dfrs.groupby('request_id')
    n_usable = grpr.agg(lambda x: np.sum(x > 0.05)).astype(int)
    n_usable.name = 'n_usable' 
    # also make a df_usable with slot columns for below
    df_usable = dfrs.unstack() > 0.05

    # sum df_metric down to one column
    metric_sum = grpr.agg(lambda x: np.sum(np.where(x > 0, x, 0)))
    metric_sum.name = 'metric_sum'

    # merge additional useful info
    dfr = df[['program_id','subprogram_name','total_requests_tonight']].join(n_usable).join(metric_sum)

    dfr['occupancy'] = dfr['total_requests_tonight']/dfr['n_usable']
    # zero out any unusable slots
    dfr.loc[dfr['n_usable'] == 0, 'occupancy'] = 0.

    # Create an empty model
    m = Model('requests')

    # set the number of threads Gurobi uses
    if 'GRB_USE_NTHREADS' in os.environ:
        m.Params.Threads = int(os.environ['GRB_USE_NTHREADS'])

    # decision variable: yes or no for each request set
    yr_dict = m.addVars(df_metric_local.index,name='Yr',vtype=GRB.BINARY)
    yr_series = pd.Series(yr_dict,name='Yr')
    dfr = dfr.join(yr_series)

    m.setObjective(np.sum(dfr['Yr'] * dfr['metric_sum'] * dfr['occupancy']), 
        GRB.MAXIMIZE)


    # slot occupancy constraint: nreqs obs divided over nusable slots
    constr_avg_slot_occupancy = m.addConstrs(
        ((np.sum(df_usable[t]*dfr['occupancy'] * dfr['Yr'])
        <= max_exps_per_slot) for t in slots), "constr_avg_slot_occupancy")

    # program balance.  To avoid key errors, only set constraints 
    # for programs that are present
    requests_needed = []
    for p in requests_allowed.keys():
        if np.sum((dfr['program_id'] == p[0]) &
                (dfr['subprogram_name'] == p[1])) > 0:
            requests_needed.append(p)

    constr_balance = m.addConstrs(
        ((np.sum(dfr.loc[(dfr['program_id'] == p[0]) & 
                (dfr['subprogram_name'] == p[1]), 'Yr'] * 
            dfr.loc[(dfr['program_id'] == p[0]) & 
                (dfr['subprogram_name'] == p[1]), 'total_requests_tonight'])
        <= requests_allowed[p]) for p in requests_needed), 
        "constr_balance")

    # Quick and dirty is okay!
    # TODO: tune this value
    m.Params.TimeLimit = 30.


    m.update()

    m.optimize()

    if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
        raise ValueError("Optimization failure")


    # now get the decision variables
    dfr['Yr_val'] = dfr['Yr'].apply(lambda x: x.getAttr('x') > 0.1)

    return dfr.loc[dfr['Yr_val'],'program_id'].index, dft


def slot_optimize(df_metric, df, requests_allowed):
    """Determine which slots to place the requests in.

    Decision variable is yes/no per request_id, slot, filter,
    with an additional decision variable on which filter to use in which slot"""

    request_sets = df_metric.index.values
    # these are fragile when columns get appended
    slots = np.unique(df_metric.columns.get_level_values(0).values)
    filter_ids = np.unique(df_metric.columns.get_level_values(1).values)
    # extra columns floating around cause problems
    filter_ids = [fid for fid in filter_ids if fid != '']

    # flatten the metric dataframe to make it easier to work with 
    df_metric_local = df_metric.copy()
    df_metric_local['request_id'] = df_metric_local.index

    # make a "tidy" dataframe with one row per (request, slot, filter)
    dft = pd.melt(df_metric_local,id_vars='request_id',
        var_name=['slot','metric_filter_id'],
        value_name='metric')
    # get n_reqs by fid
    n_reqs_cols = ['n_reqs_{}'.format(fid) for fid in filter_ids]
    n_reqs_cols.extend(['program_id','subprogram_name','total_requests_tonight'])
    dft = pd.merge(dft,df[n_reqs_cols],left_on='request_id',right_index=True)


    # Create an empty model
    m = Model('slots')

    # set the number of threads Gurobi uses
    if 'GRB_USE_NTHREADS' in os.environ:
        m.Params.Threads = int(os.environ['GRB_USE_NTHREADS'])

    yrtf_dict = m.addVars(dft.index,name='Yrtf',vtype=GRB.BINARY)
    yrtf_series = pd.Series(yrtf_dict,name='Yrtf')
    dft = dft.join(yrtf_series)

    m.setObjective(np.sum(dft['Yrtf'] * dft['metric']), 
        GRB.MAXIMIZE)


    # no more than nreqs_{fid} slots assigned per request set
    # TODO: this constructor is pretty slow
    constr_nreqs = m.addConstrs(
        ((np.sum(dft.loc[(dft['request_id'] == r) & 
                        (dft['metric_filter_id'] == f), 'Yrtf']) 
                        <= df.loc[r,'n_reqs_{}'.format(f)])
                        for f in filter_ids for r in request_sets), 
                        "constr_nreqs")

    # create resultant variables: Ytf = 1 if slot t has filter f used
    ytf = m.addVars(slots, filter_ids, vtype=GRB.BINARY)
    for t in slots:
        for f in filter_ids:
            m.addGenConstrOr(ytf[t,f], 
                dft.loc[(dft['slot'] == t) &
                        (dft['metric_filter_id'] == f), 'Yrtf'], #tolist() 
                        "orconstr_{}_{}".format(t,f))

    # now constrain ourselves to one and only one filter per slot.  
    constr_onefilter = m.addConstrs(
        (ytf.sum(t,'*') == 1 for t in slots), 'constr_onefilter')
    

    # no more than nobs per slot
    # TODO: tune this so it's closer to the average performance so we don't
    # drop too many...
    constr_nperslot = m.addConstrs(
        ((np.sum(dft.loc[dft['slot'] == t, 'Yrtf'])
        <= max_exps_per_slot) for t in slots), "constr_nperslot")

    # program balance.  To avoid key errors, only set constraints 
    # for programs that are present
    requests_needed = []
    for p in requests_allowed.keys():
        if np.sum((dft['program_id'] == p[0]) &
                (dft['subprogram_name'] == p[1])) > 0:
            requests_needed.append(p)

    constr_balance = m.addConstrs(
        ((np.sum(dft.loc[(dft['program_id'] == p[0]) & 
            (dft['subprogram_name'] == p[1]), 'Yrtf'])
        <= requests_allowed[p]) for p in requests_needed), 
        "constr_balance")

    # set a minimum metric value we'll allow, so that flagged limiting mags
    # (-99) are locked out: 1e-5 is limiting mag ~13
    # need to use generalized constraints 

    # this ought to work but gurobi can't seem to parse the sense variable
    # correctly
    #constr_min_metric = m.addConstrs((((row['Yrtf'] == 1) >> (row['metric'] >= 1.e-5)) for (_, row) in dft.iterrows()), "constr_min_metric")
    # so do the slow loop:
#    for idx, row in dft.iterrows():
#        m.addGenConstrIndicator(row['Yrtf'], True, row['metric'], 
#                GRB.GREATER_EQUAL, 1.e-5)

    # sadly above leads to infeasible models


    # Quick and dirty is okay!
    # TODO: tune this value
    m.Params.TimeLimit = 30.

    m.update()

    m.optimize()

    if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
        raise ValueError("Optimization failure")


    # now get the decision variables.  Use > a constant to avoid 
    # numerical precision issues
    dft['Yrtf_val'] = dft['Yrtf'].apply(lambda x: x.getAttr('x') > 0.1) 

    df_schedule = dft.loc[dft['Yrtf_val'],['slot','metric_filter_id', 'request_id']]

    n_iterations = 1    
    # if we don't optimize long enough, we can end up not satisfying
    # our constraints.  In that case, continue the optimization
    # TODO: verify all of them
    while df_schedule.groupby(['slot','request_id']).agg(len).max()[0] > 1:
        n_iterations += 1
        if n_iterations > 10:
            1/0
        print("> Slot optimization did not satisfy all constraints. Continuing Optimization (Iteration {})".format(n_iterations)) 
        m.update()
        m.optimize()

        # now get the decision variables
        dft['Yrtf_val'] = dft['Yrtf'].apply(lambda x: x.getAttr('x') > 0.1)
        df_schedule = dft.loc[dft['Yrtf_val'],['slot','metric_filter_id', 'request_id']]



    return df_schedule

def tsp_optimize(pairwise_distances):
    # core algorithmic code from
    # http://examples.gurobi.com/traveling-salesman-problem/

    # Callback - use lazy constraints to eliminate sub-tours 
    def subtourelim(model, where): 
        if where == GRB.callback.MIPSOL: 
            selected = [] 
            # make a list of edges selected in the solution 
            for i in range(n): 
                sol = model.cbGetSolution([model._vars[i,j] for j in range(n)]) 
                selected += [(i,j) for j in range(n) if sol[j] > 0.5] 
            # find the shortest cycle in the selected edge list 
            tour = subtour(selected) 
            if len(tour) < n: 
                # add a subtour elimination constraint 
                expr = 0 
                for i in range(len(tour)): 
                    for j in range(i+1, len(tour)): 
                        expr += model._vars[tour[i], tour[j]] 
                model.cbLazy(expr <= len(tour)-1) 

    # Given a list of edges, finds the shortest subtour 
    def subtour(edges): 
        visited = [False]*n 
        cycles = [] 
        lengths = [] 
        selected = [[] for i in range(n)] 
        for x,y in edges: 
            selected[x].append(y) 
        while True: 
            current = visited.index(False) 
            thiscycle = [current] 
            while True: 
                visited[current] = True 
                neighbors = [x for x in selected[current] if not visited[x]] 
                if len(neighbors) == 0: 
                    break 
                current = neighbors[0] 
                thiscycle.append(current) 
            cycles.append(thiscycle) 
            lengths.append(len(thiscycle)) 
            if sum(lengths) == n: 
                break 
        return cycles[lengths.index(min(lengths))] 

    assert (pairwise_distances.shape[0] == pairwise_distances.shape[1])
    n = pairwise_distances.shape[0]

    # avoid optimization failures if we only feed in a couple of points
    if n == 1:
        return [0], [READOUT_TIME.to(u.second).value]
    if n == 2:
        return [0, 1], [pairwise_distances[0,1]]
    
    m = Model() 

    # set the number of threads Gurobi uses
    if 'GRB_USE_NTHREADS' in os.environ:
        m.Params.Threads = int(os.environ['GRB_USE_NTHREADS'])
        

    # Create variables 
    vars = {} 
    for i in range(n): 
        for j in range(i+1): 
            vars[i,j] = m.addVar(obj=pairwise_distances[i,j], vtype=GRB.BINARY, name='e'+str(i)+'_'+str(j)) 
            vars[j,i] = vars[i,j] 
        m.update() 

    # Add degree-2 constraint, and forbid loops 
    for i in range(n): 
        m.addConstr(quicksum(vars[i,j] for j in range(n)) == 2) 
        vars[i,i].ub = 0 
    m.update() 
    # Optimize model 
    m._vars = vars 
    m.params.LazyConstraints = 1 
    m.optimize(subtourelim) 

    if m.Status != GRB.OPTIMAL:
        raise ValueError("Optimization failure")

    # TODO: I could cut this block if needed
    solution = m.getAttr('x', vars) 
    selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5] 
    distances = np.sum([pairwise_distances[s] for s in selected])
    distance = m.objVal
    assert len(subtour(selected)) == n

    # dictionary of connected nodes
    edges = defaultdict(list)
    for i in range(n): 
        for j in range(n): 
            if vars[i,j].getAttr('x') > 0.5:
                edges[i].append(j)

    def unwrap_tour(edges, start_node=None):
        if start_node is None:
            start_node = 0
        
        current_node = start_node 
        # arbitrary choice of direction
        next_node = edges[start_node][0]
        tour = [start_node]

        while next_node != start_node:
            tour.append(next_node)
            edge_nodes = edges[next_node]
            assert (current_node in edge_nodes)
            assert(len(edge_nodes) == 2)
            if edge_nodes[0] == current_node:
                tmp = edge_nodes[1]
            elif edge_nodes[1] == current_node:
                tmp = edge_nodes[0]
            current_node = next_node
            next_node = tmp

        return tour
            

    tour = unwrap_tour(edges)
    assert (len(tour) == n)

    return tour, distance
