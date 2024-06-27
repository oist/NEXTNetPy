#pragma once

#include <iostream>
#include <numeric>
#include <vector>
#include <chrono>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "networkx.hpp"


namespace py = pybind11;


// Wrapper to measure the time taken to run a simulation on a given network
template <typename T>
std::vector<double> measure_run_time(py::object graph,T psi, T* rho,bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS,bool VERBOSE){
    
    if (VERBOSE){
        py::print("load network");
    }

    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    if (VERBOSE)
        py::print("starting simulation....");

    std::vector<double> measured_time({});
    for (int sim = 0; sim < NB_SIMULATIONS; sim++){

        if (VERBOSE)
            py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");


        auto start = std::chrono::high_resolution_clock::now();
                
        simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

        if (VERBOSE)
            py::print("initialise sim....");
        // Initialise
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
        if (VERBOSE)
            py::print("begin infections....");      

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

        }
        if (VERBOSE)
            py::print("end epidemic");

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        const double time = duration.count() ;
        measured_time.push_back(time);
    }
    if (VERBOSE)
        py::print("leaving function");

    return measured_time;
};


// Wrapper to measure the time taken to run a simulation on a given network
template <typename T>
std::vector<double> measure_run_time_nmga(py::object graph,T psi, T* rho,bool SIR,double TMAX, int INITIAL_INFECTED, int seed, int NB_SIMULATIONS,bool VERBOSE){
    
    if (VERBOSE){
        py::print("load network");
    }

    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    if (VERBOSE)
        py::print("starting simulation....");

    std::vector<double> measured_time({});
    for (int sim = 0; sim < NB_SIMULATIONS; sim++){

        if (VERBOSE)
            py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");


        auto start = std::chrono::high_resolution_clock::now();
        
        simulate_nmga::params p;
	    p.approximation_threshold = SIZE/10;
        p.SIR = SIR;
        
        simulate_nmga simulation(network,psi,rho,p);
        // simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

        if (VERBOSE)
            py::print("initialise sim....");
        // Initialise
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
        if (VERBOSE)
            py::print("begin infections....");      

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

        }
        if (VERBOSE)
            py::print("end epidemic");

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        const double time = duration.count() ;
        measured_time.push_back(time);
    }
    if (VERBOSE)
        py::print("leaving function");

    return measured_time;
};