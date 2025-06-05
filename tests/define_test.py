import nextnet as nn
import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def test_static_SIR(show_plot = False):

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
    if show_plot:
        plt.plot(results["time"],results["infected"])
        plt.xlabel("time")
        plt.ylabel("number of infected")
        plt.show()


    if len(results["infected"]) != 2*10**5:
        print("test_static_SIR failed")
        return False
    else: return True

def test_activity_driven(show_plot = False):

    # Set up random generator
    seed = 0
    rng = nn.rng(seed)

    n = 10**3 # size of the network
    activity_rates= [1]*n
    m = 3 # number of edges added per node
    graph = nn.activity_driven_network(activity_rates, 1, m, .01, rng)

    # Define the distribution for the infection times
    psi = nn.transmission_time_exponential(10)

    # Define the simulation object
    sim = nn.simulation_temporal(graph,psi)

    # Add initial infections
    initial_infected = [(0,0)]
    sim.add_infections(initial_infected)

    # Run simulation
    opt = {"max_infected":n,"network_events":False,"epidemic_events":True}
    results = sim.run(rng,opt)

    # Plot results
    if show_plot:
        plt.plot(results["time"],results["infected"])
        plt.xlabel("time")
        plt.ylabel("number of infected")
        plt.show()

    if len(results["infected"]) != n-1:
        print("test_activity_driven failed")
        print(len(results["infected"]))
        return False
    else: return True
