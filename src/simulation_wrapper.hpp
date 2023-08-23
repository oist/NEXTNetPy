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


/** Wrapper to run a simulation in DISCRETE steps n.
* Returns a tuple:
* - the average number of newly infected Z(n), n=0,1,2....
* - the solution of the recursion equation
* - the average degree of the leaves at step n.
*/
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_leaves(py::object graph,int sim=100,int seed=1);

std::tuple<std::vector<std::vector<double>>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_zk(py::object graph,int SIZE, int sim=100, int seed=1);


// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> simulate(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1);

// Wrapper to run the average simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1, int NB_SIMULATIONS = 1, bool TRIM = false);

// Wrapper to run a simulation on a 2D LATTICE given a network and transmission distribution
// it was used to create nice GIFs to illustrate the propagation of a non-Markovian epidemic on a 2D lattice.
std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object graph,int ROWS,int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, const std::vector<double>& check_times = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0},int seed = 1);
