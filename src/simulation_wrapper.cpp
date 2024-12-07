#include <iostream>
#include <numeric>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "algorithm.h"
#include "analysis.h"
#include "temporal_network.h"
#include "network.h"
#include "networkx.hpp"
#include "NextReaction.h"
#include "random.h"
#include "REGIR.h"

#include "simulation_wrapper.hpp"
#include "tools.hpp"

std::tuple<std::vector<double>, std::vector<int>> simulate(network& network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed){
    
    rng_t engine;
    engine.seed(seed);

    const int SIZE = (int) network.nodes();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;


    simulate_next_reaction::params p;
    p.shuffle_neighbours = SHUFFLE_NEIGHBOURS;
    p.edges_concurrent = EDGES_CONCURRENT;
    p.SIR = SIR;
    simulate_next_reaction simulation(network, psi,rho,p);

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
            case epidemic_event_kind::outside_infection:
            case epidemic_event_kind::infection:
                current_number_infected++;
                break;
            case epidemic_event_kind::reset:
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

std::tuple<std::vector<double>, std::vector<int>> simulate(py::object py_nw,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed){

    rng_t engine;
    engine.seed(seed);

    networkx nw(py_nw);

    const int SIZE = (int) nw.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;


    simulate_next_reaction::params p;
    p.shuffle_neighbours = SHUFFLE_NEIGHBOURS;
    p.edges_concurrent = EDGES_CONCURRENT;
    p.SIR = SIR;
    simulate_next_reaction simulation(nw, psi,rho,p);

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
            case epidemic_event_kind::outside_infection:
            case epidemic_event_kind::infection:
                current_number_infected++;
                break;
            case epidemic_event_kind::reset:
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

std::tuple<std::vector<double>, std::vector<double>> simulate_average(py::object py_nw,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES){
    
    rng_t engine;
    engine.seed(seed);

    networkx nw(py_nw);

    const int SIZE = (int) nw.adjacencylist.size();

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

        simulate_next_reaction::params p;
        p.shuffle_neighbours = SHUFFLE_NEIGHBOURS;
        p.edges_concurrent = EDGES_CONCURRENT;
        p.SIR = SIR;
        simulate_next_reaction simulation(nw, psi,rho,p);


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
                case epidemic_event_kind::outside_infection:
                case epidemic_event_kind::infection:
                    current_number_infected++;
                    trajectory.push_back(std::make_pair(point->time,(double) 1.0/NB_SIMULATIONS) );
                    break;
                case epidemic_event_kind::reset:
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

