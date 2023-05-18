#pragma once 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "NextReaction.h"

namespace py=pybind11;

double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi);

std::tuple<std::vector<double>,double> generalised_knn(py::object graph, int moment = 1, double r = 0.0,int seed = 1);

/**
 * @brief Function that returns the empirical degree distribution of the susceptible nodes at different stages of the epidemic.
 * 
 * 
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph, transmission_time_gamma psi,const std::vector<double>& freq = {0.1,0.25,0.5,0.75,0.9},int seed = 1);


// py::list freq = py::list({0.1,0.25,0.5,0.75,0.9})

// void save_network_epidemic_state(simulate_next_reaction& simulation, std::string filename);

void save_grid(std::vector<std::vector<int>>& grid, std::string filename);