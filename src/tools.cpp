#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "tools.hpp"
#include "NextReaction.h"
#include "networkx.hpp"
#include "graph.h"
#include "nMGA.h"


namespace py=pybind11;

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>> knn_pk(int SIZE, int SIM, int seed){
    rng_t engine;
    engine.seed(seed);
    
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    transmission_time_deterministic psi(1);
    
    int n_min = 40;
    
    int global_kmax = SIZE;
    // knn_average[n] returns average knn at the n_th step of the epidemic.
    // std::vector<std::unordered_map<int,double>> knn_average(n_min);
    // std::vector<std::unordered_map<int,double>> knn_counts(n_min);
    std::vector<std::vector<double>> knn_average(n_min,std::vector<double>(SIZE,0));
    std::vector<std::vector<double>> knn_counts(n_min,std::vector<double>(SIZE,0));
    std::vector<std::vector<double>> pk_average(n_min,std::vector<double>(SIZE,0));

    for (int s = 1; s<=SIM;s++){
        py::print(s,"/",SIM,"\r",py::arg("end") = "");

        scale_free network(SIZE,engine);   
        // erdos_reyni network(SIZE,5,engine);     
        simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

        //0. obtain initial knn and pk
        std::vector<double> pk_num(SIZE,0);
        // std::vector<double> KNN_num(SIZE,0);
        // std::vector<double> KNN_den(SIZE,0);

        int local_kmax= 0;

        for (node_t node=0; node < SIZE; node++){
            const int k = network.outdegree(node);
            local_kmax = std::max(local_kmax,k);
            
            // pk
            pk_num[k] += 1.0 ;

            // knn
            for (node_t neigh : network.adjacencylist[node])
            {
                const double k_neigh = (double) network.outdegree(neigh);
                knn_average[0][k] += k_neigh;
                knn_counts[0][k] += 1.0;
            }            
        }
    
        for (int k=0;k <= local_kmax; k++){
            pk_average[0][k] += (double) pk_num[k] / (SIZE * SIM);
            // const bool ZERO = (KNN_den[k]==0);
            // knn_average[0][k] += (ZERO) ? 0 : KNN_num[k] / (KNN_den[k] * SIM);
        }
    
        // Initial infection
        std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
        const node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
        
        int n = 0;
        // int current_step = 0;
        int SIZE_LEFT = SIZE;

        // begin simulation
        while (true) {
            
            auto point = simulation.step(engine);
            if (!point)
                break;


            const node_t infected_node = point->node;
            const int current_step = (int) std::round(point->time);
            const int k = network.outdegree(infected_node);

            // update pk
            pk_num[k] --;
            SIZE_LEFT --;
            
            // keep track of knn
            for (node_t neigh : network.adjacencylist[infected_node])
            {
                const double k_neigh = (double) network.outdegree(neigh);
                knn_average[current_step][k] -= k_neigh;
                knn_counts[current_step][k] -= 1.0;
            }         

            // if generation completed, normalise pk and export knn
            if (current_step > n){
                for (int k = 0; k<=local_kmax; k++){
                    pk_average[current_step][k] += (double) pk_num[k] / (SIZE_LEFT * SIM);

                    // const bool ZERO = (KNN_den[k]==0);
                    // knn_average[n][k] += (ZERO) ? 0 : KNN_num[k] / (KNN_den[k] * SIM);
                }
                n = current_step;
            }
        
        } // epidemic ended
        // update kmax: for the average to be effective we want all arrays to be the same size, hence this min-max k-degree.
        global_kmax = std::min(local_kmax,global_kmax);

        // if the epidemic never went exponential (died out after a few steps, say, < 6) then skip that step
        // note: this part of the code is never needed in BA networks as the network is connected (i.e no disconnected components).
        if (n<=5)
            continue;
        n_min = std::min(n,n_min);

    } // simulations ended.

    // remove unused steps
    while (knn_average.size() > n_min ){
        knn_average.pop_back();
        pk_average.pop_back();
    }

    // remove unused degrees
    for (int n = 0; n < n_min; n++){
        while (knn_average[n].size() > global_kmax){
            knn_average[n].pop_back();
            pk_average[n].pop_back();
        }

        for (int k = 0; k < global_kmax; k++){
            const bool ZERO = (knn_counts[n][k]==0);
            knn_average[n][k] = ZERO ? 0 : knn_average[n][k]/knn_counts[n][k];
        }
    }

    

    //temp test
    double sum = 0;
    for (auto p : pk_average[0])
        sum += p;
    py::print("sum is ...",sum,"\n");



    return std::make_tuple(knn_average,pk_average);
}

std::vector<double> analyse_leaves(std::vector<node_t>& leaves,graph_adjacencylist& nw, simulate_next_reaction* simulate, int kmax){

    const int SIZE = (int) nw.adjacencylist.size();
    const int SIZE_LEAVES = (int) leaves.size();
    
    std::vector<double> KNN(kmax+1, 0);
    std::vector<double> KNN_counts(kmax+1,0);
    std::vector<double> PK(kmax+1,0);

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    // py::print("leaves size : ",SIZE_LEAVES,"\n");
    for (node_t node : leaves){
        
        // 1. compute the moments and the degree distribution
        int k = nw.outdegree(node);

        k1 += (double) k / SIZE_LEAVES;
        k2 += (double) pow(k,2) / SIZE_LEAVES;
        k3 += (double) pow(k,3) / SIZE_LEAVES;
        PK[k] += (double) 1 / SIZE_LEAVES ;

        // 2. compute knn. i.e. the average degree of the neighbours of a leaf with deg. k.
        // here in two parts. first compute the sum of the degree, and normalise later after the loop.
        for (node_t neigh : nw.adjacencylist[node])
        {
            // we discard the ancestor node that gave the disease to this leaf.
            if (simulate -> is_infected(neigh))
                continue;
            const double k_neigh = (double) nw.outdegree(neigh);
            KNN[k] += k_neigh;
            KNN_counts[k] += 1.0;
        }
    }

    // normalise the counts to express the correct knn
    for (int k=0; k < kmax+1; k++)
    {
        if (KNN_counts[k]!=0){
            KNN[k] /= KNN_counts[k];
        }
    }

    // 3. compute assortativity
    // sigma * r = [sum_k knn(k) * k^2 *p[k] / <k> ] - mu**2

    const double MU = k2/k1;
    const double SIGMA = k3/k1 - pow(MU,2);

    double numerator=0;
    for (int k = 0; k <= kmax; k++){
        numerator += KNN[k]*k*k*PK[k]/k1;
    }
    numerator -= pow(MU,2);

    const double r = numerator/SIGMA;
    
    return {r,k1,k2,k3};
    

}

std::vector<double> knn_depleted(graph_adjacencylist& nw, simulate_next_reaction* simulate){

    int size = (int) nw.adjacencylist.size();
    std::vector<double> knn_degree(size, 0);
    std::vector<double> nb_with_degree(size,0);

    int kmax = 0;
    for (node_t node = 0; node < size; node++){
        if (simulate && simulate -> is_infected(node))
            continue;
        int k = nw.outdegree(node);
        kmax = std::max(k,kmax);
        for (node_t neigh : nw.adjacencylist[node])
        {
            // if (simulate && simulate -> is_infected(neigh))
            //     continue;
            const double k_neigh = (double) nw.outdegree(neigh);
            knn_degree[k] += k_neigh;
            nb_with_degree[k] += 1;
        }

    }

    while ((int) knn_degree.size() > kmax + 1){
        knn_degree.pop_back();
        nb_with_degree.pop_back();
    }
    
    for (int k=0; k < kmax+1; k++)
    {
        if (nb_with_degree[k]!=0){
            knn_degree[k] =  knn_degree[k]/nb_with_degree[k];
        }
    }
    return knn_degree;

}



double assortativity_depleted(std::vector<node_t>& leaves, graph_adjacencylist& nw, simulate_next_reaction* simulate){

   
    const int SIZE = (int) nw.adjacencylist.size();
    std::vector<double> Knn = knn_depleted(nw,simulate);
    
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;

    int kmax = (int) Knn.size();
    std::vector<double> pk(SIZE,0);
    
    int SIZE_LEFT = SIZE;
    for( int node = 0; node < SIZE; node++ ){
        if(simulate && simulate -> is_infected(node)){
            SIZE_LEFT --;
            continue;
        }
        const int k = nw.outdegree(node);
        pk[k] += (double) 1 ;
        k1 += (double) k ;
        k2 += (double) k*k ;
        k3 += (double) k*k*k ;
    }

    for (int k = 0; k<= kmax; k++)
        pk[k] /= SIZE_LEFT;
    
    k1 /= SIZE_LEFT;
    k2 /= SIZE_LEFT;
    k3 /= SIZE_LEFT;
    
    while(pk.size() > kmax)
        pk.pop_back();
    
    const double MU = k2/k1;
    const double SIGMA = k3/k1 - pow(MU,2);

    double numerator=0;
    for (int k = 0; k <= kmax; k++){
        numerator += Knn[k]*k*k*pk[k]/k1;
    }
    numerator -= pow(MU,2);

    return numerator/SIGMA;




}


// measure_runtime(std::vector<int>& sizes,std::vector<double>& average_run_time, std::vector<double>& std_run_time, Factory factory,
// 					 int NB_SIM, double TMAX,int MAX_POWER, std::string filename,rng_t& engine)


// std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> benchmark(std::string network_ensemble, std::string simulation_method,transmission_time_gamma psi, transmission_time_gamma* rho,bool SIR,int BA_M, int NB_SIM, double TMAX,int MAX_POWER, bool EDGES_CONCURRENT, int INITIAL_INFECTED,int seed){


//     if (simulation_method == "NR"){

        

//     } else if (simulation_method == "NMGA"){

//     } else if (simulation_method == "REGIR"){

//     } else 



// }

// std::tuple<std::vector<int>,std::vector<double>, std::vector<double>> run_benchmark(std::string ensemble, transmission_time_gamma psi, transmission_time_gamma* rho,bool SIR,int BA_M, int NB_SIM, double TMAX,int MAX_POWER, bool EDGES_CONCURRENT, int INITIAL_INFECTED,int seed){
    
//     // declare random number generator.
//     rng_t engine(seed);

//     // import the networkx module
//     py::object nx = py::module_::import("networkx");

//     // Method that will call the appropriate networkx function
//     py::object graph_generator;

//     py::object graph;


//     const bool SHUFFLE_NEIGHBOURS = false;
//     const double MEAN = 5;
//     const double VARIANCE = 3;

//     transmission_time_gamma psi(MEAN, VARIANCE);

//     std::vector<double> average_run_time({});
//     std::vector<double> std_run_time({});
//     std::vector<int> network_size({});

//     for (int power = 7; power < 21; power ++){
//         const int SIZE = pow(2,power);

//         py::print("N = ", power - 6,"/", 15,"\r",py::arg("end") = "");

//         network_size.push_back(SIZE);


//         // Determine the network ensemble and generate a new network
//         if (ensemble=="ER"){
//                 graph_generator = nx.attr("fast_gnp_random_graph");
//                 graph = graph_generator(SIZE,(double) 3/SIZE);

//         }
//         else if (ensemble=="BA"){
//             graph_generator = nx.attr("barabasi_albert_graph");
//             graph = graph_generator(SIZE,BA_M);
//         }
//         else if (ensemble == "WS"){
//             graph_generator = nx.attr("watts_strogatz_graph");
//             graph = graph_generator(SIZE,2,0.15);
//             graph = nx.attr("convert_node_labels_to_integers")(graph);
//         }
//         else {
//             throw std::logic_error("unrecognised ensemble name\n The possible names are ER, BA, WS");
//         }

//         networkx network(graph);

//         // To sample random nodes uniformly
//         std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

//         double mean = 0;
//         double mean_square = 0;
// 		for (int i = 0; i < NB_SIM; i++)
// 		{

// 			// Create simulation environment
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

//     return std::make_tuple(network_size, average_run_time, std_run_time);
    
// }


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

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph,const std::vector<double>& freq,int seed){

    std::vector<double> fraction;
    for (auto fi : freq){
        fraction.push_back(fi);
        // py::print("frac ",fi,"\n");
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


    transmission_time_deterministic psi(1);

    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;

    simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT);

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
            // py::print("found ",f,"\n");
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


std::vector<double> euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi){
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

    return {GROWTH_RATE,MU} ;
}