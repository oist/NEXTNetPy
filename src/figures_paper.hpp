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

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> simulate_on_powerlaw(int SIZE,double EXPONENT,transmission_time_gamma psi, transmission_time_gamma* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed,int NB_SIMULATIONS, bool TRIM);

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> simulate_on_lognormal(int SIZE,double mean, double variance,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1,int NB_SIMULATIONS=1, bool TRIM=false);