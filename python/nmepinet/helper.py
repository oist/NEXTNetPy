import nmepinet as nmen
import numpy as np
import time
from tqdm import tqdm
import networkx as nx

def measure_time(N, nb_simulations,a=0.1,b=1):
    activity_rates = [a] * N
    inactivation_rate =b 
    eta = 1
    m = 3
    beta = 1
    mu = 0.1
    tmax = 5000
    seed = 1
    initial_infected = 10
    t0 = 0

    times = []

    for s in tqdm(range(nb_simulations)):
        start_time = time.time()
        
        time_traj,infected = nmen.simulate_on_activity_driven(
                activity_rates,
                inactivation_rate,
                eta,
                m,
                beta,
                mu,
                TMAX=tmax,
                seed=seed+s,
                initial_infected=initial_infected,
                t0=t0,
                nb_simulations=1,
                SIR=True
                )
        
        end_time = time.time()
        times.append(end_time - start_time)

    # Calculate mean and standard deviation
    avg_time = np.mean(times)
    std_dev_time = np.std(times)

    return avg_time, std_dev_time

def measure_time_static(N, nb_simulations,a=0.1,b=1):
    m = 3
    beta = 1
    psi= nmen.transmission_time_exponential(beta)
    mu = 0.1
    rho= nmen.transmission_time_exponential(mu)
    tmax = 5000
    seed = 1
    initial_infected = 10

    graph = nx.fast_gnp_random_graph(N,10/N)

    times = []

    for s in tqdm(range(nb_simulations)):
        start_time = time.time()
        
        time_traj,infected = nmen.simulate(
                graph,
                psi,
                rho,
                TMAX=tmax,
                seed=seed+s,
                initial_infected=initial_infected,
                SIR=True
                )
        
        end_time = time.time()
        times.append(end_time - start_time)

    # Calculate mean and standard deviation
    avg_time = np.mean(times)
    std_dev_time = np.std(times)

    return avg_time, std_dev_time


def simulate_on_activity_driven2(graph, psi,rho=None,SIR=True,initial_infected =1, nb_simulations=1,seed=1,verbose=False,trim=True,tmax=1000):

    # graph = nmen.activity_driven_graph(activity_rates,eta,m,recovery_rate,rng)
    og = graph
    
    time_trajectory=[]
    infected_trajectory=[]

    for s in range(nb_simulations):
        if verbose:
            print(str(s)+'/'+str(nb_simulations),end='\r')

        graph = og 
    
        time,infected = nmen.simulate_on_temporal(
            graph,
            psi,
            rho,
            SIR=SIR,
            initial_infected=initial_infected,
            seed=seed+s,
            verbose=verbose,
            network_size=graph.nodes(),
            TMAX=tmax,
            trim=True)
                
    
        time_trajectory +=time
        infected_trajectory += infected

    # Convert to numpy arrays
    time_trajectory = np.array(time_trajectory)
    infected_trajectory = np.array(infected_trajectory)

    # Get the sorted indices based on time_trajectory
    sorted_indices = np.argsort(time_trajectory)

    # Sort both arrays using the sorted indices
    time_trajectory = time_trajectory[sorted_indices]
    if  nb_simulations > 1 and trim:
        return time_trajectory[::nb_simulations], np.cumsum(infected_trajectory[sorted_indices])[::nb_simulations]/nb_simulations
    else:
        return time_trajectory, np.cumsum(infected_trajectory[sorted_indices])/nb_simulations




def simulate_SIRX(graph,kappa0,kappa, psi,rho=None,SIR=True,initial_infected =1, nb_simulations=1,seed=1,verbose=False,trim=True,tmax=1000):


    time_trajectory=[]
    infected_trajectory=[]
    if verbose:
        print("begin simulation:")
    for s in range(nb_simulations):
        if verbose:
            print(str(s)+'/'+str(nb_simulations),end='\r')

        temporal_graph =nmen.dynamic_sirx_network(graph,kappa0,kappa)
        time,infected = nmen.simulate_on_temporal(
            temporal_graph,
            psi,
            rho,
            SIR=SIR,
            initial_infected=initial_infected,
            seed=seed+s,
            verbose=verbose,
            network_size=graph.nodes(),
            TMAX=tmax,
            trim=True)
                
    
        time_trajectory +=time
        infected_trajectory += infected

    # Convert to numpy arrays
    time_trajectory = np.array(time_trajectory)
    infected_trajectory = np.array(infected_trajectory)

    # Get the sorted indices based on time_trajectory
    sorted_indices = np.argsort(time_trajectory)

    # Sort both arrays using the sorted indices
    time_trajectory = time_trajectory[sorted_indices]
    if  nb_simulations > 1 and trim:
        return time_trajectory[::nb_simulations], np.cumsum(infected_trajectory[sorted_indices])[::nb_simulations]/nb_simulations
    else:
        return time_trajectory, np.cumsum(infected_trajectory[sorted_indices])/nb_simulations


def simulate_on_empirical(path,psi,rho=None,finite_duration=False,dt=1,SIR=True,initial_infected =1, nb_simulations=1,seed=1,verbose=False,trim=True,tmax=1000):

    number_of_nodes = 0  # Start with a very low value
    
    with open(path, 'r') as file:  # Replace with your file path
        for line in file:
            first_column = int(line.split()[0])  # Split the line and 
            number_of_nodes = max(first_column,number_of_nodes)
    
    time_trajectory=[]
    infected_trajectory=[]
    if verbose:
        print("found ",number_of_nodes,"nodes")
        print("begin simulation:")
    for s in range(nb_simulations):
        if verbose:
            print(str(s)+'/'+str(nb_simulations),end='\r')

        temporal_graph = nmen.temporal_empirical_graph(path,finite_duration,dt)
    
        time,infected = nmen.simulate_on_temporal(
            temporal_graph,
            psi,
            rho,
            SIR=SIR,
            initial_infected=initial_infected,
            seed=seed+s,
            verbose=verbose,
            network_size=number_of_nodes,
            TMAX=tmax,
            trim=True)
                
    
        time_trajectory +=time
        infected_trajectory += infected

    # Convert to numpy arrays
    time_trajectory = np.array(time_trajectory)
    infected_trajectory = np.array(infected_trajectory)

    # Get the sorted indices based on time_trajectory
    sorted_indices = np.argsort(time_trajectory)

    # Sort both arrays using the sorted indices
    time_trajectory = time_trajectory[sorted_indices]
    if  nb_simulations > 1 and trim:
        return time_trajectory[::nb_simulations], np.cumsum(infected_trajectory[sorted_indices])[::nb_simulations]/nb_simulations
    else:
        return time_trajectory, np.cumsum(infected_trajectory[sorted_indices])/nb_simulations
