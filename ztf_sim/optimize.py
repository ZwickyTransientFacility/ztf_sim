"""Implementation of core scheduling algorithms using Gurobi."""

import logging
import os
from collections import defaultdict
from gurobipy import *
import numpy as np
import shelve
import astropy.units as u
import pandas as pd
from collections import defaultdict
from .constants import TIME_BLOCK_SIZE, EXPOSURE_TIME, READOUT_TIME, FILTER_CHANGE_TIME
from .constants import PROGRAM_NAME_TO_ID, slew_time
from .Fields import Fields
from .utils import RA_to_HA, HA_to_RA, block_index

max_exps_per_slot = np.ceil((TIME_BLOCK_SIZE / 
                (EXPOSURE_TIME + READOUT_TIME)).to(
                u.dimensionless_unscaled).value).astype(int)

logger = logging.getLogger(__name__)

def night_optimize(df_metric, df, requests_allowed, time_limit=30*u.second,
        block_use=defaultdict(float)):
    """Assign request sets to time slots using integer linear programming.

    Formulates and solves a binary ILP that assigns each request set to one
    or more 30-minute time slots, maximising a
    weighted observability metric subject to time-capacity, program-balance,
    and single-filter-per-slot constraints.

    Parameters
    ----------
    df_metric : pandas.DataFrame
        Multi-level column DataFrame indexed by ``request_set_id``. Columns
        are ``(slot, filter_id)`` pairs; values are the per-(request, slot,
        filter) priority metric.
    df : pandas.DataFrame
        Request pool DataFrame indexed by ``request_set_id``. Required
        columns: ``program_id``, ``subprogram_name``,
        ``total_requests_tonight``, ``exposure_time``, ``dec``, and
        ``n_reqs_{fid}`` for each filter ID.
    requests_allowed : dict
        Mapping ``(program_id, subprogram_name) -> int`` giving the maximum
        (or exact, for MSIP) number of slot-observations allowed tonight.
    time_limit : astropy.units.Quantity, optional
        Wall-clock time budget for the Gurobi solver. Default is 30 s.
    block_use : collections.defaultdict, optional
        Mapping ``slot_index -> float`` giving the fraction of each slot
        already consumed by timed programmes. Default is all zeros.

    Returns
    -------
    request_sets_tonight : pandas.Index
        IDs of request sets selected by the optimiser.
    df_slots : pandas.DataFrame
        Tidy assignment table with columns ``'slot'``,
        ``'metric_filter_id'``, and ``'request_id'`` for every
        selected (request, slot, filter) triple.
    dft : pandas.DataFrame
        Full ``(request, slot, filter)`` DataFrame including the binary
        decision variable values (``'Yrtf_val'``).

    Raises
    ------
    ValueError
        If the optimiser cannot satisfy constraints after 10 iterations.

    Notes
    -----
    Decision variables:

    * ``Yr[r]`` — binary; request set *r* is observed tonight.
    * ``Yrtf[r,t,f]`` — binary; *r* assigned to slot *t* with filter *f*.
    * ``Ytf[t,f]`` — binary; filter *f* used in slot *t*.
    * ``Ydfds[s]`` — binary; filter change between slot *s* and *s+1*.

    Objective: maximise the weighted metric sum minus a filter-change
    penalty. MSIP programme balance is enforced with an equality constraint
    (relaxed to ≤ if infeasible); other programmes use an upper bound.
    """

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
    n_reqs_cols.extend(['program_id','subprogram_name',
        'total_requests_tonight','exposure_time','dec'])
    dft = pd.merge(dft,df[n_reqs_cols],left_on='request_id',right_index=True)

    # don't need the dec column anymore
    dft = dft.drop('dec',axis=1)

    # calculate number of slots required per request set
    
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

    # sum df_metric down to one column
    metric_sum = grpr.agg(lambda x: np.sum(np.where(x > 0, x, 0)))
    metric_sum.name = 'metric_sum'

    # merge additional useful info
    dfr = df[['program_id','subprogram_name','total_requests_tonight']].join(n_usable).join(metric_sum)

    # determine which request sets have enough usable slots
    dfr['observable_tonight'] = dfr['total_requests_tonight'] <= dfr['n_usable']
 
    # restrict to only the observable requests
    dfr = dfr.loc[dfr['observable_tonight'],:]
    dft = pd.merge(dft,dfr[['n_usable','observable_tonight']],
            left_on='request_id',right_index=True)
    dft = dft.loc[dft['observable_tonight'],:]
    request_sets = dfr.index.values
    df_metric = df_metric.loc[dfr.index]

    # Create an empty model
    m = Model('slots')

    # set the number of threads Gurobi uses
    if 'GRB_USE_NTHREADS' in os.environ:
        m.Params.Threads = int(os.environ['GRB_USE_NTHREADS'])

    # decision variable: yes or no for each request set
    yr_dict = m.addVars(df_metric_local.index,name='Yr',vtype=GRB.BINARY)
    yr_series = pd.Series(yr_dict,name='Yr')
    dfr = dfr.join(yr_series)


    yrtf_dict = m.addVars(dft.index,name='Yrtf',vtype=GRB.BINARY)
    yrtf_series = pd.Series(yrtf_dict,name='Yrtf')
    dft = dft.join(yrtf_series)

    # putting the Yrtf vars in the dft series is not always convenient;
    # provide a lookup table
    rtf_to_idx = {}
    idx_to_rtf = dft[['request_id','slot','metric_filter_id']].to_dict(orient='index')
    for k,v in idx_to_rtf.items():
        rtf_to_idx[(v['request_id'],v['slot'],v['metric_filter_id'])] = k
    


    # create resultant variables: Yr = 1 if request r is observed in at least
    # one slot
    for r in request_sets:
        m.addGenConstrOr(yr_dict[r], dft.loc[dft['request_id'] == r, 'Yrtf'],
                "orconstr_{}".format(r))

    # nreqs_{fid} slots assigned per request set if it is observed
    # this constructor is pretty slow
    constr_nreqs = m.addConstrs(
        ((np.sum(dft.loc[(dft['request_id'] == r) & 
                        (dft['metric_filter_id'] == f), 'Yrtf']) 
                        == (df.loc[r,'n_reqs_{}'.format(f)] * dfr.loc[r,'Yr']))
                        for f in filter_ids for r in request_sets), 
                        "constr_nreqs")

    # minimum slot separation per filter constraint
    MIN_SLOT_SEPARATION = 2
    # TODO: generalize this beyond just ZUDS
    wZUDSt = ( (dft['subprogram_name'] == 'high_cadence') 
              )

    # only add these parameters if there is a program to space 
    if np.sum(wZUDSt):
        space_obs = True
    else:
        space_obs = False

    if space_obs:
        wZUDSg = wZUDSt & (dft['metric_filter_id'] == 1)
        wZUDSr = wZUDSt & (dft['metric_filter_id'] == 2)
        wZUDSi = wZUDSt & (dft['metric_filter_id'] == 3)
        wrZUDS =  ((dft['subprogram_name'] == 'high_cadence'))  
        ZUDS_request_sets = dfr.loc[wrZUDS].index.tolist()
        filter_ids_to_limit = [1,2]

        # slot separation

        # create resultant variables: 1 if both slot_a and slot_b are true
        yrttf = m.addVars(ZUDS_request_sets,slots[:-1],slots[1:],
                filter_ids_to_limit, vtype=GRB.BINARY)
        constraint_dict = {}
        for r in ZUDS_request_sets:
            for t in slots[:-1]:
                for t2 in slots[1:]:
                    if t2 < t: 
                        # avoid duplicate entries
                        continue
                    #TODO: this should be removed if we're not doing hard constraints
    #                dt = t2 - t
    #                if dt >= MIN_SLOT_SEPARATION:
    #                    continue
                    for f in filter_ids_to_limit:
                        m.addGenConstrAnd(yrttf[r,t,t2,f], 
                            [yrtf_dict[rtf_to_idx[(r,t,f)]], 
                             yrtf_dict[rtf_to_idx[(r,t2,f)]]],
                            "slotdiff_and_{}_{}_{}_{}".format(r,t,t2,f))

        dtdict = defaultdict(list)
        for t in slots[:-1]:
            for t2 in slots[1:]:
                if t2 <= t: 
                    # avoid duplicate entries
                    continue
                dt = t2 - t
                dtdict[dt].append((t,t2))

        # create delta-t resultant variables: OR constraint for all pairwise slots
        yrdtf = m.addVars(ZUDS_request_sets,dtdict.keys(),
                filter_ids_to_limit, vtype=GRB.BINARY)
        for r in ZUDS_request_sets:
            for dt in dtdict.keys():
                for f in filter_ids_to_limit:
                        # loop over items in dtdict
                    m.addGenConstrOr(yrdtf[r,dt,f], 
                       [yrttf[r,t,t2,f] for (t,t2) in dtdict[dt]],
                        "slot_dt_indicator_{}_{}_{}".format(r,dt,f))

    # THIS WORKS to set hard slot separation constraints
    ##    # can set a hard constraint here by requiring all the low-separation pairs
    ##    # to be zero
    #    constr_min_slotsep = m.addConstrs(
    #        (yrdtf[r,dt,f] == 0 for r in ZUDS_request_sets for dt in dtdict.keys() if dt <= (MIN_SLOT_SEPARATION-1) for f in filter_ids_to_limit), 'constr_min_slot_sep')
    ##    # or use in the objective function
    
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

    # create filter change resultant variable: Ydfds = 1 if
    # filter changes between slot s and s+1
    ydfds = m.addVars(slots[:-1], vtype=GRB.BINARY)
    # use indicator constraints to set the value
    for i,t in enumerate(slots[:-1]):
        for f in filter_ids:
            m.addGenConstrIndicator(ydfds[t], False,
                ytf[slots[i],f] - ytf[slots[i+1], f], GRB.EQUAL, 0,
                        "filt_change_indicator_{}_{}".format(t,f))
    

    # total exposure time constraint 
    constr_nperslot = m.addConstrs(
        ((np.sum(dft.loc[dft['slot'] == t, 'Yrtf'] * 
            (dft.loc[dft['slot'] == t, 'exposure_time'] + 
                READOUT_TIME.to(u.second).value))
            <= (TIME_BLOCK_SIZE.to(u.second).value * (1. - block_use[t]))) 
            for t in slots), "constr_nperslot")

    # program balance.  To avoid key errors, only set constraints 
    # for programs that are present
    msip_requests_needed = []
    msip_requests_possible = {} 
    requests_needed = []
    for p in requests_allowed.keys():
        if p[0] == PROGRAM_NAME_TO_ID['MSIP']:
            wmsipp = (dfr['program_id'] == p[0]) & (dfr['subprogram_name'] == p[1])
            n_available = np.sum(dfr.loc[wmsipp,'total_requests_tonight'])
            if n_available > 0:
                # to demand exact equality we need to know how many 
                # requests we have
                # set to the minimum of allowed or available
                if n_available <=  requests_allowed[p]:
                    # MSIP requests only come in pairs, so we need an even
                    # number.
                    # TODO: generalize this.
                    if n_available % 2 != 0:
                        n_available -= 1
                    msip_requests_possible[p] = n_available

                else:
                    if requests_allowed[p] % 2 != 0:
                        requests_allowed[p] += 1
                    msip_requests_possible[p] = requests_allowed[p]
                msip_requests_needed.append(p)
        else: 
            if np.sum((dft['program_id'] == p[0]) &
                    (dft['subprogram_name'] == p[1])) > 0:
                requests_needed.append(p)

    # demand exact equality for MSIP
    constr_msip_balance = m.addConstrs(
        ((np.sum(dft.loc[(dft['program_id'] == p[0]) & 
            (dft['subprogram_name'] == p[1]), 'Yrtf'])
        == msip_requests_possible[p]) for p in msip_requests_needed), 
        "constr_msip_balance")

    constr_balance = m.addConstrs(
        ((np.sum(dft.loc[(dft['program_id'] == p[0]) & 
            (dft['subprogram_name'] == p[1]), 'Yrtf'])
        <= requests_allowed[p]) for p in requests_needed), 
        "constr_balance")

    m.update()

    # np.heaviside returns a TypeError so make our own
    def heaviside(x, x0=0):
        # scalars only
        # < and > are not implimented for Gurobi Linexps, so have to do 
        # some unusual control flow here with ==, <=, >=
        if x == 0:
            return x0
        elif x <= 0:
            return 0
        else:
            return 1

    # scale by number of standard exposures so long exposures aren't
    # penalized
    if not space_obs:
        m.setObjective(
            np.sum(dft['Yrtf'] * dft['metric'] * 
            dft['exposure_time']/EXPOSURE_TIME.to(u.second).value) 
            - ydfds.sum() * (FILTER_CHANGE_TIME / (EXPOSURE_TIME +
                READOUT_TIME) * 2.5).value,
# TODO: fails under Gurobi 9.0+ ("Constraint has no bool value")
#            - np.sum(
#                [heaviside((requests_allowed[p] - np.sum(
#                    dft.loc[(dft['program_id'] == p[0]) &
#                    (dft['subprogram_name'] == p[1]), 'Yrtf'].values
#                    )))*2.5 
#                    for p in requests_needed]),
            GRB.MAXIMIZE)
    else:
        def slot_scale(dt):
            return dt/24.

        m.setObjective(
            np.sum(dft['Yrtf'] * dft['metric'] * 
            dft['exposure_time']/EXPOSURE_TIME.to(u.second).value) 
            - ydfds.sum() * (FILTER_CHANGE_TIME / (EXPOSURE_TIME +
                READOUT_TIME) * 2.5).value
            + np.sum(yrdtf[r,dt,f]*slot_scale(dt) for r in ZUDS_request_sets for dt in dtdict.keys() if dt >= MIN_SLOT_SEPARATION for f in filter_ids_to_limit), 
# TODO: fails under Gurobi 9.0+ ("Constraint has no bool value")
#            - np.sum(
#                [heaviside((requests_allowed[p] - np.sum(
#                    dft.loc[(dft['program_id'] == p[0]) &
#                    (dft['subprogram_name'] == p[1]), 'Yrtf'].values
#                    )))*2.5 
#                    for p in requests_needed]),
            GRB.MAXIMIZE)

    # Quick and dirty is okay!
    m.Params.TimeLimit = time_limit.to(u.second).value

    m.update()

    m.optimize()

    # if we have optimization problems, MSIP exact equality is likely the problem.
    # relax the constraint and retry
    

    if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
        logger.warning(f"Gurobi failed to optimize! Code {m.Status}")
        logger.info("Relaxing MSIP exact constraint.")
        m.setAttr("Sense", [c for c in constr_msip_balance.values()],
                           ["<" for c in constr_msip_balance.values()]) 
        m.optimize()
        if (m.Status != GRB.OPTIMAL) and (m.Status != GRB.TIME_LIMIT):
            logger.error(f"Gurobi optimization with relaxed MSIP constraints failed! Code {m.Status}")


    # now get the decision variables.  Use > a constant to avoid 
    # numerical precision issues
    try:
        dft['Yrtf_val'] = dft['Yrtf'].apply(lambda x: x.getAttr('x') > 0.1) 
    except AttributeError:
        logger.error("Optimization reached time limit but didn't find solutions")
        logger.info("Relaxing MSIP exact constraint.")
        m.setAttr("Sense", [c for c in constr_msip_balance.values()],
                           ["<" for c in constr_msip_balance.values()]) 
        m.optimize()
        try:
            dft['Yrtf_val'] = dft['Yrtf'].apply(lambda x: x.getAttr('x') > 0.1) 
        except AttributeError:
            logger.error(f"Gurobi optimization with relaxed MSIP constraints failed!")


    df_schedule = dft.loc[dft['Yrtf_val'],['slot','metric_filter_id', 'request_id']]

    n_iterations = 1    
    # if we don't optimize long enough, we can end up not satisfying
    # our constraints.  In that case, continue the optimization
    while df_schedule.groupby(['slot','request_id']).agg(len).max()[0] > 1:
        n_iterations += 1
        if n_iterations > 10:
            raise ValueError('Optimization failed to satisfy constraints')
        print("> Slot optimization did not satisfy all constraints. Continuing Optimization (Iteration {})".format(n_iterations)) 
        m.update()
        m.optimize()

        # now get the decision variables
        dft['Yrtf_val'] = dft['Yrtf'].apply(lambda x: x.getAttr('x') > 0.1)
        df_schedule = dft.loc[dft['Yrtf_val'],['slot','metric_filter_id', 'request_id']]

    # get the request set decision variables
    dfr['Yr_val'] = dfr['Yr'].apply(lambda x: x.getAttr('x') > 0.1)

    # this doesn't work in the objective function but is a useful check
    def num_filter_changes(ytf):

        n_changes = 0
        for i, slot in enumerate(slots[:-1]):
            for fid in filter_ids:
                if ytf[(slot,fid)].getAttr('x') == 1:
                    if not (ytf[(slots[i+1], fid)].getAttr('x') == 1):
                        n_changes+=1
        return n_changes

    print(f'Number of filter changes: {num_filter_changes(ytf)}')

    return dfr.loc[dfr['Yr_val'],'program_id'].index, df_schedule, dft

def tsp_optimize(pairwise_distances):
    """Solve the Travelling Salesman Problem for a set of telescope pointings.

    Finds the minimum-cost Hamiltonian cycle through *n* nodes using a
    Gurobi MIP formulation with lazy subtour-elimination constraints.

    Parameters
    ----------
    pairwise_distances : numpy.ndarray, shape (n, n)
        Symmetric matrix of pairwise travel costs (e.g. slew times in
        seconds). Diagonal should be zero.

    Returns
    -------
    tour : list of int
        Node visit order as a list of 0-based indices of length *n*.
    distance : float
        Total tour cost.

    Raises
    ------
    ValueError
        If the Gurobi optimisation fails to find an optimal solution.

    Notes
    -----
    Degree-2 constraint: each node is incident to exactly two edges.
    Subtour elimination is enforced lazily via a
    ``GRB.callback.MIPSOL`` callback (``subtourelim``): whenever the
    solver finds a solution containing a proper sub-cycle of length < n,
    a cut is added requiring at most ``|subtour| - 1`` edges within that
    subset. Special cases: n = 1 returns ``([0], [READOUT_TIME])``;
    n = 2 returns the single edge without calling Gurobi.
    """
    # core algorithmic code from
    # http://examples.gurobi.com/traveling-salesman-problem/

    # Callback - use lazy constraints to eliminate sub-tours
    def subtourelim(model, where):
        """Gurobi lazy-constraint callback that adds subtour elimination cuts."""
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
        """Find the shortest cycle in a list of directed edges."""
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
        """Convert an edge adjacency dict into an ordered tour list."""
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


def check_limits_and_solve_TSP(unique_fields, time_start, time_end):
    """Check telescope pointing limits and solve the TSP for a standalone queue.

    Intended for time-windowed queues such as Einstein Probe observations.
    Validates each field against the P48 pointing limits at both the start
    and end of the window, builds a pairwise slew-time matrix, prepends a
    fake starting node at the CALSTOW position (HA = 0, Dec = −48°), solves
    the TSP, then removes the fake node from the result.

    Parameters
    ----------
    unique_fields : list of int or array-like of int
        ZTF field IDs to sequence.
    time_start : astropy.time.Time
        Start of the observation window.
    time_end : astropy.time.Time
        End of the observation window.

    Returns
    -------
    queue_order : numpy.ndarray of int
        Field IDs in optimised visit order (fake start node removed).
    """

    # this is somewhat expensive to construct each time
    fld = Fields()

    # check for telescope limits
    for field in unique_fields:
        row = fld.fields.loc[int(field)]
        ra = row['ra']
        dec = row['dec']

        # TODO: this should be logged rather than printed
        if dec < -22:
            for t, limit in zip([time_start, time_end], ['start','end']):
                ha = RA_to_HA(ra*u.degree,t).to(u.degree).value
                if ha <= -17.6:
                    print(f'{field} violates limit A at {t} [{limit}]')
                if (dec < -45) & (ha >= -17.6):
                    print(f'{field} violates limit A at {t} [{limit}]')
                if (dec < -46) & (np.abs(ha) >= 3.):
                    print(f'{field} violates limits B {t} [{limit}]')


    # solve the TSP.  
    df = fld.fields.loc[list(unique_fields)]
    fld.compute_blocks(time_start)
    queue_slot = block_index(time_start)[0]

    az = fld.block_az[queue_slot]
    az.name = 'azimuth'
    df = df.join(az, on='field_id')
    # Need to use correct HA=0
    df_blockstart = pd.DataFrame({'ra':HA_to_RA(0,
        time_start).to(u.degree).value,
        'dec':-48.,'azimuth':180.},index=[0])
    df_fakestart = pd.concat([df_blockstart,df],sort=True)

    # below copied from QueueManager.py, consider a refactor
    # compute pairwise slew times by axis for all pointings
    slews_by_axis = {}
    def coord_to_slewtime(coord, axis=None):
        c1, c2 = np.meshgrid(coord, coord)
        dangle = np.abs(c1 - c2)
        angle = np.where(dangle < (360. - dangle), dangle, 360. - dangle)
        return slew_time(axis, angle * u.deg)

    slews_by_axis['dome'] = coord_to_slewtime(
        df_fakestart['azimuth'], axis='dome')
    slews_by_axis['dec'] = coord_to_slewtime(
        df_fakestart['dec'], axis='dec')
    slews_by_axis['ra'] = coord_to_slewtime(
        df_fakestart['ra'], axis='ha')

    maxradec = np.maximum(slews_by_axis['ra'], slews_by_axis['dec'])
    maxslews = np.maximum(slews_by_axis['dome'], maxradec)
    # impose a penalty on zero-length slews (which by construction
    # in this mode are from different programs)
    wnoslew = maxslews == 0
    maxslews[wnoslew] = READOUT_TIME * 10.
    overhead_time = np.maximum(maxslews, READOUT_TIME)

    tsp_order, tsp_overhead_time = tsp_optimize(overhead_time.value)

    # remove the fake starting point.  tsp_optimize always starts with
    # the first observation in df, which by construction is our fake point,
    # so we can simply cut it off.
    tsp_order = tsp_order[1:]

    # tsp_order is 0-indexed from overhead time, so I need to
    # reconstruct the request_id
    queue_order = df_fakestart.index.values[tsp_order]

    return queue_order
