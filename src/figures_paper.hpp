#pragma once

#include <iostream>
#include <numeric>
#include <vector>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "networkx.hpp"



namespace py = pybind11;

// figure 3a : power law network with disassortativity
template <typename Transmission_time>
std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> simulate_on_powerlaw(int SIZE,double EXPONENT,Transmission_time psi, Transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed,int NB_SIMULATIONS, bool TRIM){

    rng_t engine;
    engine.seed(seed);

    py::print("creating network...\r",py::arg("end") = "");
    std::vector<int> degrees = powerlaw_degree_list(EXPONENT,SIZE,engine);
    config_model network(degrees,engine);
    
    // add_correlation(-0.1,network,engine);
    
    py::print("computing moments...");

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double r = assortativity(network);
    py::print("FINAL assortativity: ",r);

    for (node_t node = 0; node < SIZE; node++){   
        const int k = network.outdegree(node);
        k1 += (double) k/SIZE;
        k2 += (double) k*k/SIZE;
        k3 += (double) k*k*k/SIZE;
    }

    std::vector<double> params = {k1,k2,k3,r};

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::vector<std::pair<double,double>> trajectory;
    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    for (int sim = 0; sim < NB_SIMULATIONS; sim++){

        py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");

        simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

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

    py::print("sorting...");

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

    return std::make_tuple(time_trajectory, infected_trajectory,params);
};


template <typename Transmission_time>
std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> simulate_on_lognormal(int SIZE,double mean, double variance,Transmission_time psi, Transmission_time* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1,int NB_SIMULATIONS=1, bool TRIM=false,bool VERBOSE=false){
    rng_t engine;
    engine.seed(seed);
    if (VERBOSE)
        py::print("creating network...\r",py::arg("end") = "");
    std::vector<int> degrees = lognormal_degree_list(mean,variance,SIZE,engine);

    config_model network(degrees,engine);
    if (VERBOSE)
        py::print("computing moments...");

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double r = assortativity(network);

    for (node_t node = 0; node < SIZE; node++){   
        const int k = network.outdegree(node);
        k1 += (double) k/SIZE;
        k2 += (double) k*k/SIZE;
        k3 += (double) k*k*k/SIZE;
    }

    std::vector<double> params = {k1,k2,k3,r};

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

    if (VERBOSE)
        py::print("sorting...");

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

    return std::make_tuple(time_trajectory, infected_trajectory,params);
};