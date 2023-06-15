#pragma once

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "analysis.h"
#include "algorithm.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "networkx.hpp"



namespace py = pybind11;

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> depletion_discrete(py::object graph,int sim=100,int seed=1);


std::tuple<std::vector<double>,std::vector<double>> simulation_discrete(py::object graph,int sim=100,int seed=1);

// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> run_simulation(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1);

// Wrapper to run the average simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= false,int INITIAL_INFECTED=1, int seed = 1, int NB_SIMULATIONS = 1, bool TRIM = false);

// Wrapper to run a simulation on a 2D LATTICE given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object graph,int ROWS,int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, const std::vector<double>& check_times = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0},int seed = 1);


// // Factory Template
// template<typename T>
// class Factory {
// public:
//     static T* createSimulation(std::string network_ensemble,int SIZE, transmission_time_gamma psi, transmission_time_gamma* rho,bool SHUFFLE_NEIGHBOURS, rng_t engine) {
        
//         if (network_ensemble == "ER"){
//             erdos_reyni network(SIZE, R0, engine);
//             return new T(network,psi,*rho,SHUFFLE_NEIGHBOURS);
//         } else if (network_ensemble == "BA"){
//             scale_free network(SIZE, engine);
//             return new T(network,psi,*rho,SHUFFLE_NEIGHBOURS);
//         }
//     }
// };