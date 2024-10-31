#pragma once

#include <iostream>
#include <numeric>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"

#include "NextReaction.h"
// #include "dynamic_graph.h"
#include "nMGA.h"
#include "networkx.hpp"
#include "REGIR.h"


namespace py = pybind11;


std::tuple<std::vector<double>, std::vector<int>> simulate(graph& network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed);
std::tuple<std::vector<double>, std::vector<int>> simulate(py::object graph,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed);

std::tuple<std::vector<double>, std::vector<double>> simulate_on_temporal(dynamic_network& network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT, int seed,int initial_infected,int size,bool trim,bool verbose,double t0);


// Wrapper to run a simulate given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> simulate_average(py::object graph,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES);
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(graph& network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES);


// Wrapper to run a simulation given a network and transmission distribution
template <typename T>
std::tuple<std::vector<double>, std::vector<int>> simulate_traj_regir(py::object graph,T psi, T* rho, bool SIR,double TMAX,int INITIAL_INFECTED, int seed){

    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    simulate_regir::params p;
	p.approximation_threshold = INITIAL_INFECTED;
    p.SIR = SIR;

	simulate_regir simulate(network, psi, rho, p);

    std::vector<double> time_trajectory({});
    std::vector<int> infected_trajectory({});

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    std::unordered_set<node_t> selected_nodes;
    for (node_t i = 0; i < INITIAL_INFECTED; i++)
    {
        const node_t random_node = uniform_node_distribution(engine);
        // sample without replacement:
        if (selected_nodes.find(random_node) != selected_nodes.end()){
            i--;
            continue;
        }
        simulate.add_infections({ std::make_pair(random_node, 0.0)});
        selected_nodes.insert(random_node);
    }
    
    int current_number_infected = 0;
                        
    while (true) {
        auto point = simulate.step(engine);
        if (!point )
            break;

        switch (point->kind) {
            case event_kind::outside_infection:
            case event_kind::infection:
                current_number_infected++;
                break;
            case event_kind::reset:
                current_number_infected--;
                break;
            default:
                break;
        }
        if (current_number_infected<=0 || point->time > TMAX)
            break;
        time_trajectory.push_back(point->time);
        infected_trajectory.push_back(current_number_infected);
    }
    return std::make_tuple(time_trajectory, infected_trajectory);
}



/** Wrapper to run a simulation in DISCRETE steps n.
* Returns a tuple:
* - the average number of newly infected Z(n), n=0,1,2....
* - the solution of the recursion equation
* - the average degree of the leaves at step n.
*/
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_leaves(py::object graph,int sim=100,int seed=1);

std::tuple<std::vector<std::vector<double>>,std::vector<int>> simulation_discrete(py::object graph,int sim=1, int seed=1,bool VERBOSE=false);


// Wrapper to run a simulation on a 2D LATTICE given a network and transmission distribution
// it was used to create nice GIFs to illustrate the propagation of a non-Markovian epidemic on a 2D lattice.
std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object graph,int ROWS,int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, const std::vector<double>& check_times = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0},int seed = 1);

// Wrapper to run a simulation given a network and transmission distribution
template <typename T>
std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<double>,std::vector<double>,std::vector<int>> simulation_per_degree_class(py::object graph,T psi, T* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES){
    
    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;


    /* Aim: create a map to get all different degrees.
    * A vector is inappropriate for degree distribution with large tails
    * because too many elements will be zero, requiring too much memory.
    */

    std::vector<int> unique_degrees({});

    double k_mean = 0; 
    for (node_t node = 0; node < SIZE; node++){
        const int k0 = network.outdegree(node);
        k_mean += (double) k0;
        auto it = std::lower_bound(unique_degrees.begin(), unique_degrees.end(), k0);
        if (it == unique_degrees.end() || *it != k0)
            unique_degrees.insert(it, k0);
    }

    const int klen = (int) unique_degrees.size();
   // Create an unordered map to store positions
    std::unordered_map<int, int> pos;
    for (int i = 0; i < klen; ++i) {
        const int k = unique_degrees[i];
        pos[k] = i;
    }

    std::vector<std::vector<std::pair<double,double>>> trajectory(klen,std::vector<std::pair<double,double>>({}));
    std::vector<std::pair<double,double>> k_mean_traj;

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    for (int sim = 0; sim < NB_SIMULATIONS; sim++){

        if (VERBOSE)
            py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");

        simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);



        if (ALL_NODES){
            const node_t initial_node = sim;
            simulation.add_infections({ std::make_pair(initial_node, 0)});
        } else {

            std::unordered_set<node_t> selected_nodes;
            for (node_t i = 0; i < INITIAL_INFECTED; i++)
            {
                const node_t random_node = uniform_node_distribution(engine);
                // sample without replacement:
                if (selected_nodes.find(random_node) != selected_nodes.end()){
                    i--;
                    continue;
                }
                simulation.add_infections({ std::make_pair(random_node, 0)});
                selected_nodes.insert(random_node);
            }
        }
        
        int current_number_infected = 0;
        double k1 = k_mean;
        double s1 = (double) SIZE;

        while (true) {
            auto point = simulation.step(engine);
            if (!point )
                break;
            
            const int k0 = network.outdegree(point->node);
            s1 --;
            k1 -= (double) k0;
            k_mean_traj.push_back(std::make_pair(point->time,k1/s1));

            const int i0 = pos[k0];

            switch (point->kind) {
                case event_kind::outside_infection:
                case event_kind::infection:
                    current_number_infected++;
                    trajectory[i0].push_back(std::make_pair(point->time,(double) 1.0/NB_SIMULATIONS) );
                    break;
                case event_kind::reset:
                    current_number_infected--;
                    trajectory[i0].push_back(std::make_pair(point->time,(double) -1.0/NB_SIMULATIONS) );
                    break;
                default:
                    break;
            }
            if (current_number_infected<=0 || point->time > TMAX)
                break;

            
        }
    }


    // All simulations ended, sort the vector (first element by default)
    std::vector<std::vector<double>> time_trajectory(klen,std::vector<double>({}));
    std::vector<std::vector<double>> infected_trajectory(klen,std::vector<double>({}));

    for (int i = 0; i < klen; i++){
        std::sort(trajectory[i].begin(),trajectory[i].end(),[](const auto& x, const auto& y){return x.first < y.first;});

        // trim the trajectory and split the vectors
        int counter = TRIM ? NB_SIMULATIONS : 1;
        double infected = 0;
        for (index_t j = 0; j < trajectory[i].size(); j++){
            infected += trajectory[i][j].second;
            if (j % counter== 0){
                time_trajectory[i].push_back(trajectory[i][j].first);
                infected_trajectory[i].push_back(infected);
            }
        }
    }

    // average degree of healthy nodes
    std::vector<double> k_traj({});
    std::vector<double> k_traj_time({});
    std::sort(k_mean_traj.begin(),k_mean_traj.end(),[](const auto& x, const auto& y){return x.first < y.first;});
    int counter = TRIM ? NB_SIMULATIONS : 1;
    double infected = 0;
    for (index_t j = 0; j < k_mean_traj.size(); j++){
        if (j % counter== 0){
            k_traj.push_back(k_mean_traj[j].second);
            k_traj_time.push_back(k_mean_traj[j].first);
        }
    }

    return std::make_tuple(time_trajectory, infected_trajectory,k_traj_time,k_traj,unique_degrees);
};
