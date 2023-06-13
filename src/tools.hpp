#pragma once 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "simulation_wrapper.hpp"

namespace py=pybind11;

std::vector<double> knn_depleted(graph_adjacencylist& nw, simulate_next_reaction* simulate=nullptr);
double assortativity_depleted(graph_adjacencylist& network, simulate_next_reaction* simulate = nullptr);


/**
 * @brief measure the average time taken to simulate an epidemic on a network ensemble.
 * 
 * 
 */
std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> run_benchmark(std::string ensemble, std::string method, transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR= true, int MAX_POWER = 21, int NB_SIM = 500, double TMAX = 1000, bool EDGES_CONCURRENT= false,int INITIAL_INFECTED = 1, int seed = 1);


std::vector<double> euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi);

std::tuple<std::vector<double>,double> generalised_knn(py::object graph, int moment = 1, double r = 0.0,int seed = 1);

/**
 * @brief Function that returns the empirical degree distribution of the susceptible nodes at different stages of the epidemic.
 * 
 * 
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph,const std::vector<double>& freq = {0.1,0.25,0.5,0.75,0.9},int seed = 1);


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
std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> run_benchmark_next_reaction(std::string network, transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR= true, double TMAX = 1000, bool EDGES_CONCURRENT= false,int INITIAL_INFECTED = 1, int seed = 1);


typedef std::function<simulation_algorithm* (rng_t& engine, int n)> simulation_factory_t;

// void measure_runtime(rng_t& engine, simulation_factory_t factory, int SMAX, double TMAX, std::string filename);

// template<typename Factory>
// void measure_runtime(std::vector<int>& sizes,std::vector<double>& average_run_time, std::vector<double>& std_run_time, Factory factory,
// 					 int NB_SIM, double TMAX,int MAX_POWER, std::string filename,rng_t& engine)
// {
	

//     for (int power = 7; power < MAX_POWER; power ++){
//         const int SIZE = pow(2,power);
// 		network_size.push_back(SIZE);

//         py::print("N = ",SIZE,", # " power - 6,"/", 15,"\r",py::arg("end") = "");

        
//         double mean = 0;
//         double mean_square = 0;
// 		for (int i = 0; i < NB_SIM; i++)
// 		{
// 			// Create simulation environment

// 			Product* productA = Factory<ConcreteProductA>::createProduct();
//     		productA->operation();

// 			auto s = factory(engine, N);
// 			simulation_algorithm& simulation = *s.simulator;

// 			// Start measuring performance
// 			auto start = std::chrono::high_resolution_clock::now();
	

// 			// Initial number of infected
//             for (node_t i = 0; i < INITIAL_INFECTED; i++)
//             {
//                 const node_t random_node = uniform_node_distribution(engine);
//                 // sample without replacement:
//                 if (simulation.is_infected(random_node)){
//                     i--;
//                     continue;
//                 }
//                 simulation.add_infections({ std::make_pair(random_node, 0)});
//             }
			
// 			// Run simulation, collect transmission times
// 			while (true) {

// 				auto point = simulation.step(engine);
// 				if (!point || (point -> time > TMAX))
// 					break;
// 			}
		
// 			auto stop = std::chrono::high_resolution_clock::now();
			
// 			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
// 			const double time = duration.count() ;
            
//             mean += time / NB_SIM ;
//             mean_square += time * time / NB_SIM ;
// 		}
//         const double std_dev = std::sqrt( mean_square - mean * mean ) ;
//         average_run_time.push_back(mean);
//         std_run_time.push_back(std_dev);
// 	}

// }