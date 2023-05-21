#include <iostream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "tools.hpp"
#include "NextReaction.h"
#include "networkx.hpp"
#include "nMGA.h"


namespace py=pybind11;

std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> run_benchmark_next_reaction(network_ensemble network, transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR= true,int BA_M=1, double TMAX = 1000, bool EDGES_CONCURRENT= false, int seed = 1){
    
    // declare random number generator.
    rng_t engine(seed);

    // import the networkx module
    py::object nx = py::module_::import("networkx");

    // Method that will call the appropriate networkx function
    py::object graph_generator;

    py::object graph;

    const bool SHUFFLE_NEIGHBOURS = false;
    const bool EDGES_CONCURRENT = true;

    std::vector<double> average_run_time({});
    std::vector<double> std_run_time({});
    std::vector<int> network_size({});



    for (int power = 7; power < 21; power ++){
        const int SIZE = pow(2,power);
        network_size.push_back(SIZE);

        // Determine the network ensemble and generate a new network
        switch(network){
            case network_ensemble::erdos_reyni: { 
                graph_generator = nx.attr("fast_gnp_random_graph");
                graph = graph_generator(SIZE,(double) 3/SIZE);
                break;
            }
            case network_ensemble::barabasi_albert: {
                graph_generator = nx.attr("barabasi_albert_graph");
                graph = graph_generator(SIZE,BA_M);
                break;
            }
            case network_ensemble::watts_strogatz: {
                graph_generator = nx.attr("watts_strogatz_graph");
                graph = graph_generator(SIZE,2,0.15);
                graph = nx.attr("convert_node_labels_to_integers")(graph);
                break;
            }
            default: {
                throw std::logic_error("unrecognised ensemble name");
                break;
            }
        }

        networkx network(graph);

        double mean = 0;
        double mean_square = 0;

        for (int simulation = 0; simulation < 500; simulation++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            
            simulate_next_reaction simulate(network,psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

            while (true) {
                auto point = simulate.step(engine);
                if (!point)
                    break;     
            }

            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            
            mean += (double) duration.count() / 500 ;
            mean_square += (double) duration.count() * duration.count() / 500 ;
            
        }

        //Compute standard deviation
        const double std_dev = std::sqrt( mean_square - mean * mean ) ;
        average_run_time.push_back(mean);
        std_run_time.push_back(std_dev);

    }

    return std::make_tuple(network_size, average_run_time, std_run_time);
    
}



void save_grid(std::vector<std::vector<int>>& grid, std::string filename){

    // Open the file stream
    std::ofstream out_file(filename);

    // Loop over the outer vector and write each inner vector to the file
    for (const auto& row : grid) {
        for (const auto& element : row) {
            out_file << element << " ";
        }
        out_file << std::endl;
    }

    // Close the file stream
    out_file.close();
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph, transmission_time_gamma psi,const std::vector<double>& freq,int seed){

    // Convert the python list into a c++ vector
    // std::vector<int> fraction(py::len(freq)); // create a C++ vector with the same size as the Python list
    // for (auto i = 0; i < py::len(freq); i++)
    //     fraction[i] = py::cast<int>(freq[i]);


    std::vector<double> fraction;
    for (auto fi : freq){
        fraction.push_back(fi);
        py::print("frac ",fi,"\n");
    }
    // vector that contains prob. deg. distr. of the sus. nodes when a fraction f of the network is infected.
    std::vector<std::vector<double>> Prob_degree_K_depleted;

    // vector that contains the evolution of the first four moments of the susceptible nodes through out the epidemic.
    std::vector<std::vector<double>> evolving_moments;
    
    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    int k1 = 0;
    int k2 = 0;
    int k3 = 0;
    int k4 = 0;

    // non-normalised current deg. distr. of the sus. nodes
    std::vector<int> degree_counts(SIZE,0);
    int k_max = 0;

    for (node_t node = 0; node < SIZE; node++ ){
        const int k = network.outdegree(node);
        k_max = std::max(k,k_max);
        degree_counts[k] += 1;
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3) ;
        k4 += pow(k,4) ;
    }

    //Delete the useless entries
    while ((int) degree_counts.size() > k_max + 1)
        degree_counts.pop_back();


    // Compute the empirical distribution before the epidemic starts.
    std::vector<double> initial_pk;
    for (int i = 0; i < k_max+1; i++){
            initial_pk.push_back((double) degree_counts[i] / SIZE);
    }
    Prob_degree_K_depleted.push_back(initial_pk);

    std::vector<double> k1_traj({(double) k1 / SIZE});
    std::vector<double> k2_traj({(double) k2 / SIZE});
    std::vector<double> k3_traj({(double) k3 / SIZE});
    std::vector<double> k4_traj({(double) k4 / SIZE});



    simulate_next_reaction simulation(network, psi);

    //Infect the first individual by choosing a node at random
    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
    const node_t random_node = uniform_node_distribution(engine);
    simulation.add_infections({ std::make_pair(random_node, 0)});

    int SIZE_LEFT = SIZE;

    // The epidemic begins
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;

        SIZE_LEFT --;
        const double f = std::ceil((double) (SIZE-SIZE_LEFT)/SIZE * 100) / 100;

        const node_t infected_node = point->node;
        const int k = network.outdegree(infected_node);


        degree_counts[k] -= 1;
        k1 -= k ;
        k2 -= pow(k,2) ;
        k3 -= pow(k,3) ;
        k4 -= pow(k,4) ;

        k1_traj.push_back((double) k1 / SIZE_LEFT);
        k2_traj.push_back((double) k2 / SIZE_LEFT);
        k3_traj.push_back((double) k3 / SIZE_LEFT);
        k4_traj.push_back((double) k4 / SIZE_LEFT);
        
        // Use std::find to check if the value is in the vector
        auto it = std::find(fraction.begin(), fraction.end(), f);
        if (it != fraction.end()) {
            py::print("found ",f,"\n");
            fraction.erase(it);

            std::vector<double> pk;
            for (int i = 0; i < k_max+1; i++){
                    pk.push_back((double) degree_counts[i] / SIZE_LEFT);
            }

            Prob_degree_K_depleted.push_back(pk);
        }

    }//  <-- Epidemic ended

    evolving_moments.push_back(k1_traj);
    evolving_moments.push_back(k2_traj);
    evolving_moments.push_back(k3_traj);
    evolving_moments.push_back(k4_traj);

    return std::make_tuple(Prob_degree_K_depleted, evolving_moments);
}


double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi){
// double euler_lotka_growth_rate(py::object graph,double mean, double variance){
    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    double r = assortativity(network);

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;

    for (node_t node = 0; node < SIZE; node++ ){
        double k = (double) network.outdegree(node);
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3) ;
    }
    const double MU = ( (1-r) * (k2/k1 - 1 ) + r *(k3/k2 - 1) );
    const double SHAPE = psi.mean * psi.mean / psi.variance;
    const double SCALE = psi.variance / psi.mean;
    // const double SHAPE = mean * mean / variance;
    // const double SCALE = variance / variance;

    //Only valid for Gamma distribution 
    const double GROWTH_RATE = 1/SCALE * ( pow(MU,1/SHAPE) - 1 );

    return GROWTH_RATE;
}