#pragma once

#include <iostream>
#include <numeric>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "networkx.hpp"



namespace py = pybind11;

// Wrapper to run a simulation given a network and transmission distribution
template <typename T>
std::tuple<std::vector<double>, std::vector<int>> simulate(py::object graph,T psi, T* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed){

    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

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
        simulation.add_infections({ std::make_pair(random_node, 0)});
        selected_nodes.insert(random_node);
    }
    
    int current_number_infected = 0;
                        
    while (true) {
        auto point = simulation.step(engine);
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



// Wrapper to run a simulation given a network and transmission distribution
template <typename T>
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,T psi, T* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES){
    
    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::vector<std::pair<double,double>> trajectory;
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
                            
        while (true) {
            auto point = simulation.step(engine);
            if (!point )
                break;

            switch (point->kind) {
                case event_kind::outside_infection:
                case event_kind::infection:
                    current_number_infected++;
                    trajectory.push_back(std::make_pair(point->time,(double) 1.0/NB_SIMULATIONS) );
                    break;
                case event_kind::reset:
                    current_number_infected--;
                    trajectory.push_back(std::make_pair(point->time,(double) -1.0/NB_SIMULATIONS) );
                    break;
                default:
                    break;
            }
            if (current_number_infected<=0 || point->time > TMAX)
                break;

            
        }
    }


    // All simulations ended, sort the vector (first element by default)
    std::sort(trajectory.begin(),trajectory.end(),[](const auto& x, const auto& y){return x.first < y.first;});

    std::vector<double> time_trajectory({});
    std::vector<double> infected_trajectory({});
        
    // trim the trajectory and split the vectors
    int counter = TRIM ? NB_SIMULATIONS : 1;
    double infected = 0;
    for (index_t i = 0; i < trajectory.size(); i++){
        infected += trajectory[i].second;
        if (i % counter== 0){
            time_trajectory.push_back(trajectory[i].first);
            infected_trajectory.push_back(infected);
        }
    }
    
    return std::make_tuple(time_trajectory, infected_trajectory);
};

/** Wrapper to run a simulation in DISCRETE steps n.
* Returns a tuple:
* - the average number of newly infected Z(n), n=0,1,2....
* - the solution of the recursion equation
* - the average degree of the leaves at step n.
*/
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_leaves(py::object graph,int sim=100,int seed=1);

std::tuple<std::vector<std::vector<double>>,std::vector<double>> simulation_discrete(py::object graph,int sim=1, int seed=1,bool VERBOSE=false);


// Wrapper to run a simulation on a 2D LATTICE given a network and transmission distribution
// it was used to create nice GIFs to illustrate the propagation of a non-Markovian epidemic on a 2D lattice.
std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object graph,int ROWS,int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, const std::vector<double>& check_times = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0},int seed = 1);
