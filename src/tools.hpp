#pragma once 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"

namespace py=pybind11;

double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi);



/**
 * @brief Function that returns the empirical degree distribution of the susceptible nodes at different stages of the epidemic.
 * 
 * 
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph, transmission_time_gamma psi,const std::vector<double>& freq = {0.1,0.25,0.5,0.75,0.9},int seed = 1);


// py::list freq = py::list({0.1,0.25,0.5,0.75,0.9})

// void save_network_epidemic_state(simulate_next_reaction& simulation, std::string filename);

void save_grid(std::vector<std::vector<int>>& grid, std::string filename);



/**simulate_next_reaction simulation(graph, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
 * @brief Type of network ensemble
 */
enum class network_ensemble : unsigned int {
    erdos_reyni = 0,
	barabasi_albert = 1,
	watts_strogatz = 2,
	barabasi_albert_5 = 3
};



/**
 * @brief Type of simulation
 */
enum class simulation_method : unsigned int {
    next_reaction = 0,
	nmga = 1,
	regir = 2
};


/**
 * @brief measure the average time taken to simulate an epidemic on a network ensemble.
 * 
 * 
 */
std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> run_benchmark_next_reaction(network_ensemble network, transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR= true, double TMAX = 1000, bool EDGES_CONCURRENT= false, int seed = 1);

