import nextnet as nn
import networkx as nx
import matplotlib.pyplot as plt

n = 10**5 # size of the network
m = 1 # number of edges added per node
graph = nx.barabasi_albert_graph(n,m)

# Define the distribution for the infection times (Gamma distributed)
MEAN_INFECTION = 5
VARIANCE_INFECTION = 1
psi = nn.transmission_time_gamma(MEAN_INFECTION,VARIANCE_INFECTION)

# Define the distribution for the recovery times
MEAN_RECOVERY = 14
VARIANCE_RECOVERY= 3
rho = nn.transmission_time_lognormal(MEAN_RECOVERY,VARIANCE_RECOVERY)

# To simulate a SIR epidemic with one initial infected individual
sim = nn.simulation(nn.networkx(graph),psi,rho,SIR=True)
sim.add_infections([(0,0)])
results = sim.run(engine=nn.rng(1))

plt.plot(results["time"], results["infected"])
plt.savefig("out/test.pdf")
