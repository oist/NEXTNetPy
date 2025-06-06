# Simulations on Static Networks

## Simulation object

The `simulation` class wraps the `simulate_next_reaction` C++ simulator for static networks.

## Constructor

```python
simulation(network, infection_time, recovery_time=None, SIR=True)
```

* **network** (`network`): your static network object
* **infection\_time** (`transmission_time`): distribution of infection delays
* **recovery\_time** (`transmission_time`, optional): distribution of recovery delays (omit for pure SI)
* **SIR** (`bool`, default `True`): choose SIR (`True`) or SIS (`False`) dynamics

Keeps the Python objects alive while the simulator exists (`keep_alive<1,2>`, etc.).

#### Methods

* **`add_infections(infections_list)`**
  Add initial “seeds” or outside infections.

  * `infections_list` (`List[Tuple[int, float]]`): sequence of `(node_index, infection_time)`.

* **`run(engine, options=None) → dict`**
  Advance the simulation until one of:

  1. simulation time exceeds `options["time"]`
  2. cumulative infections exceed `options["total_infected"]`
  3. concurrent infections exceed `options["max_infected"]`
  4. no infected remain

  **Arguments:**

  * `engine` (`rng`): a random‐number generator instance (see example)
  * `options` (`dict`, optional): see “Options keys” above

  **Returns:**
  A dict containing:

  * `"time"`: event times
  * `"infected"`: current infected counts
  * `"recovered"`: current recovered counts
  * `"data"`: raw event tuples `(time, node, source_node, event_type)` where `event_type` is `0` for infection and `1` for recovery/reset.

### Example

```python
import nextnet as nn
import numpy as np

# Set random generator
seed = 0
engine = nn.rng(seed)

# Build a network (e.g. Erdős–Rényi)
nw = nn.erdos_renyi(1000, 4,engine)

# Define transmission/recovery distributions
rate_lambda = 0.5
rate_mu = 0.1
infect_dist = transmission_time.exponential(rate_lambda)
recovery_dist = transmission_time.exponential(rate_mu)

# Instantiate the simulator (SIR)
sim = nn.simulation(nw, infect_dist, recovery_dist, SIR=True)

# 4. seed initial infections at time 0
initial = [(i, 0.0) for i in np.random.choice(nw.size(), 5, replace=False)]
sim.add_infections(initial)

# 5. run until time 100 or 1000 total infections

result = sim.run(engine, {"time": 100.0, "total_infected": 1000})

# 6. inspect results
times = result["time"]
infected = result["infected"]
recovered = result["recovered"]
```
