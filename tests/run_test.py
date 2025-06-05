import nextnet as nn
import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# import matplotlib.pyplot as plt

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
sim = nn.simulation(graph,psi,rho,SIR=True)

# Add initial infections
t0=0
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