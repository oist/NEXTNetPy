#pragma once

#include <iostream>
#include <numeric>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"

#include "NextReaction.h"
// #include "temporal_network.h"
#include "nMGA.h"
#include "networkx.hpp"
#include "REGIR.h"


namespace py = pybind11;


std::tuple<std::vector<double>, std::vector<int>> simulate(network& network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed);
std::tuple<std::vector<double>, std::vector<int>> simulate(py::object network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed);

//simulate_on_activity_average disables because it depends on average_trajectories which
//is only meant for unit tests and doesn't work with the bundled boost
//std::tuple<std::vector<double>, std::vector<double>> simulate_on_activity_average(std::vector<double>& activity_rates, double inactivation_rate,double eta, int m, double beta, double mu,double TMAX,int seed,int initial_infected,double t0,int nb_simulations,bool SIR);

// Wrapper to run a simulate given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> simulate_average(py::object network,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool ALL_NODES);


// Wrapper to run a simulation given a network and transmission distribution
template <typename T>
std::tuple<std::vector<double>, std::vector<int>> simulate_traj_regir(py::object py_network,T psi, T* rho, bool SIR,double TMAX,int INITIAL_INFECTED, int seed){

    rng_t engine;
    engine.seed(seed);

    networkx nw(py_network);

    const int SIZE = (int) nw.adjacencylist.size();

    simulate_regir::params p;
	p.approximation_threshold = INITIAL_INFECTED;
    p.SIR = SIR;

	simulate_regir simulate(nw, psi, rho, p);

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



/** Wrapper to run a simulation in DISCRETE steps n.
* Returns a tuple:
* - the average number of newly infected Z(n), n=0,1,2....
* - the solution of the recursion equation
* - the average degree of the leaves at step n.
*/
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_leaves(py::object network,int sim=100,int seed=1);

std::tuple<std::vector<std::vector<double>>,std::vector<int>> simulation_discrete(py::object network,int sim=1, int seed=1,bool VERBOSE=false);


// Wrapper to run a simulation on a 2D LATTICE given a network and transmission distribution
// it was used to create nice GIFs to illustrate the propagation of a non-Markovian epidemic on a 2D lattice.
std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object network,int ROWS,int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, const std::vector<double>& check_times = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0},int seed = 1);
