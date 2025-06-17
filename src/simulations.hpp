#pragma once

#include <iostream>
#include <numeric>
#include "nextnet/algorithm.h"
#include "nextnet/random.h"
#include "nextnet/NextReaction.h"

#include "networkx.hpp"


// Wrapper to run a simulate given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> simulate_average(py::object graph,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM,bool VERBOSE, bool INITIAL_WITH_BIAS);
