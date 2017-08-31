
# Queues

The scheduler supports several queue methods.

Simple list queue manager: directly specify a sequence of observations, 
with optional expiration times.

Greedy queue: out of a request set of observations, take the next best one, 
either using a metric derived from volumetric survey speed 
or (to be implemented) a user-supplied metric.

Gurobi queue: a slot-based queue that optimizes observations through the 
entire night.
