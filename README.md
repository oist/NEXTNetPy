# NEXT-NetPy

*NEXT-NetPy* is a python module can simulate non-Markovian epidemics on arbitrary networks. NEXT-NetPy is based on a pybind11 wrapper around the C++ simulator https://github.com/oist/NEXTNet

## Set up

### From Release
Download the latest release and in a new terminal window run:

```bash
pip install nextnet-0.2.0.tar.gz
```
### From Source

In a new terminal window run:
```bash
#Download the repository
git clone git@github.com:oist/NEXTNetPy.git
cd NEXTNetPy
# Download the required dependencies (boost, pybind11)
git submodule update --init --recursive
# install the python package
pip install .
```

## Getting started

Here is an example of how to run a SIR simulation on a Barabasi-Albert random graph using networkx.

```python
import nextnet as nn
import networkx as nx
import matplotlib.pyplot as plt

n = 10**5 # size of the network
m = 1 # number of edges added per node
graph_nx = nx.barabasi_albert_graph(n,m)

# convert network to a NextNet object.
graph = nn.networkx(graph_nx)


# Define the distribution for the infection times (Gamma distributed)
MEAN_INFECTION = 5
VARIANCE_INFECTION = 1
psi = nn.transmission_time_gamma(MEAN_INFECTION,VARIANCE_INFECTION)

# Define the distribution for the recovery times
MEAN_RECOVERY = 14
VARIANCE_RECOVERY= 3
rho = nn.transmission_time_lognormal(MEAN_RECOVERY,VARIANCE_RECOVERY)

# Define the simulation object
sim = nn.simulate(graph,psi,rho,SIR=True)

# Add initial infections (node 0 infected at time t=0)
initial_infected = [(0,0)]
sim.add_infections(initial_infected)

# Set up random generator
seed = 0
rng = nn.rng(seed)

# Run simulation
results = sim.run(engine=rng)

# Plot results
plt.plot(results["time"],results["infected"])
plt.xlabel("time")
plt.ylabel("number of infected")
plt.show()
```

<p align="center">
  <img src="images/example_SIR.png" alt="Alt text for the image" width="300"/>
</p>


## Documentation
https://samuelcure.github.io/NextNetPy-Doc/
