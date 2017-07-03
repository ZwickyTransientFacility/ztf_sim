

from gurobipy import *
import numpy as np
import shelve
import astropy.units as u
import pandas as pd
from collections import defaultdict
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
