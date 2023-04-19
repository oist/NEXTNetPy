#pragma once 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"

namespace py=pybind11;

double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi);



/**
 * @brief Function that returns the empirical degree distribution of the susceptible nodes at different stages of the epidemic.
 * 
 * 
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph, transmission_time_gamma psi,const std::vector<double>& freq = {0.1,0.25,0.5,0.75,0.9},int seed = 1);


// py::list freq = py::list({0.1,0.25,0.5,0.75,0.9})