#include <iostream>
#include <numeric>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"
#include "networkx.hpp"
#include "simulation_wrapper.hpp"
#include "tools.hpp"
#include "graph.h"

std::tuple<std::vector<std::vector<double>>,std::vector<double>> simulation_discrete(py::object graph,int sim, int seed,bool VERBOSE){
    rng_t engine;
    engine.seed(seed);
    
    if (VERBOSE)
        py::print("reading network...");

    // convert the python graph into our c++ network object
    networkx network(graph);
    // erdos_reyni network(10000,4,engine);
    const int SIZE = (int) network.adjacencylist.size();

    if (VERBOSE)
        py::print("computing moments...");
    // Define the degree distribution and the moments
    std::vector<double> pk(SIZE,0);
    int kmax = 0;
    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double r = assortativity(network);
    
    for (node_t node = 0; node < SIZE; node++){   
        const int k = network.outdegree(node);
        pk[k] += (double) 1 / SIZE;
        kmax = std::max(kmax,k);
        k1 += (double) k/SIZE;
        k2 += (double) k*k/SIZE;
        k3 += (double) k*k*k/SIZE;
    }
    while (pk.size()>kmax+1)
        pk.pop_back();

    // Parameters of the simulation method
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    // the epidemic is in discrete steps
    transmission_time_deterministic psi(1);
    
    // Define Z[n,k] as the number of nodes that get infected at step n and have degree k
    // (Even for tree-networks, n rarely goes above 15 in epidemics)
    int n_min = 30;
    std::vector<std::vector<double>> zn_average(n_min,std::vector<double>(kmax+1,0));
    
    // mu_n is same as k2_traj/k1_traj -1 but instead of being a fct of the # of infected, it is a fct of n
    std::vector<double> mu_n(n_min,0);
    mu_n[0]=sim*(k2/k1-1);

    // use to normalise over the number of simulations. 
    // Necessary because each simulation may not reach the same number of generations: some epidemics finish in 4 steps while others may finish in 8 steps.
    // Therefore, need to normalise accordingly.
    std::vector<double> step_counter(n_min,0);
    
    if (VERBOSE)
        py::print("begin epidemic...");
    // Used to sample the initial infected node uniformly
    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    // Beginning of the simulations
    for (int s = 0; s<sim;s++){
        py::print(s,"/",sim,"\r",py::arg("end") = "");
        
        simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
       
        std::unordered_set<node_t> selected_nodes;
        for (node_t i = 0; i < 1; i++)
        {
            const node_t random_node = uniform_node_distribution(engine);
            // sample without replacement:
            if (selected_nodes.find(random_node) != selected_nodes.end()){
                i--;
                continue;
            }
            simulation.add_infections({ std::make_pair(random_node, 0)});
            selected_nodes.insert(random_node);
        }
             
        int n = 0;
        int recorded_step = 0;

        int infected = 0;

        double k1_t = k1;
        double k2_t = k2;

        if (VERBOSE)
            py::print("current step",recorded_step);
        
        int nb_infected_in_curr_gen = 0;
        while (true) {

            auto point = simulation.step(engine);
            if (!point)
                break;

            infected ++;
            const node_t infected_node = point->node;
            const int k = network.outdegree(infected_node);


            n = (int) std::round(point->time);

            zn_average[n][k] ++ ;

            // if we reach the next generation of infected update counter
            if (n>recorded_step){
                if (VERBOSE)
                    py::print("rechead new point",n);
                mu_n[n] += (k2_t / k1_t) - 1 ;
                step_counter[recorded_step]++;
                recorded_step = n;
            }

            // This is updated AFTER mu because mu_n represents the state before the infected of the nth generation contribute.
            k1_t -= (double) k/SIZE;
            k2_t -= (double) k*k/SIZE;

        }
        if (VERBOSE)
            py::print("epidemic ended");
        // Epidemic ended, but still need to update the last step
        // Note that we use n instead of recorded_step in case there is only one infected in gen n.
        step_counter[n]++;
    }

    //delete unused steps
    while (step_counter.back()==0 ){
        step_counter.pop_back();
        zn_average.pop_back();
        mu_n.pop_back();
    }

    //normalise zn:
    for (int n = 0; n<step_counter.size(); n++){
        const double norm = (double) step_counter[n];
        for( int k =0;k<= kmax;k++){
            zn_average[n][k] /= norm;
        }
    }

    //normalise mu_n:
    for (int n = 0; n<step_counter.size(); n++){
        if (mu_n[n]==0)
            continue;
        mu_n[n]/=step_counter[n];
    }
    if (VERBOSE)
        py::print("success");
    return std::make_tuple(zn_average,mu_n);
}




// rng_t engine;
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> simulation_discrete_leaves(py::object graph,int sim, int seed){
    rng_t engine;
    engine.seed(seed);
    networkx network(graph);
    
    const int SIZE = (int) network.adjacencylist.size();
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    transmission_time_deterministic psi(1);
    
    int n_min = 50;
    std::vector<double> zn_average(n_min,0);
    std::vector<double> average_leaf_degree(n_min,0);
    
    for (int s = 0; s<sim;s++){
        py::print(s,"/",sim,"\r",py::arg("end") = "");
        
        simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
        std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
       
        const node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
        
        int current_step = 0;
        int current_infected = 0;
        int n = 0;
        double zn = 1;
        int contributions = 0;
        int norm = 0;
        // begin simulation
        while (true) {
            
            auto point = simulation.step(engine);
            if (!point)
                break;

            const node_t infected_node = point->node;
            n = (int) std::round(point->time);

            const int k = network.outdegree(infected_node);
            average_leaf_degree[n] += (double) k / sim ;

            zn_average[n] += (double) 1 / sim ;
           
        }


        if (n<=5)
            continue;
        n_min = std::min(n,n_min);
    }
    while (zn_average.size() > n_min ){
        zn_average.pop_back();
        average_leaf_degree.pop_back();
    }

    for (int i =0;i<average_leaf_degree.size();i++){
        average_leaf_degree[i] /= zn_average[i];
    }

    std::vector<double> recursion(n_min,0);
    recursion[0]=1;
    recursion[1]=average_leaf_degree[0];

    for (int i =2;i<average_leaf_degree.size();i++){
        recursion[i]= recursion[i-1]*(average_leaf_degree[i-1]-1); 
    }
    return std::make_tuple(zn_average,recursion,average_leaf_degree);
}



std::tuple<std::vector<double>, std::vector<int>> run_simulation_lattice(py::object graph,int ROWS, int COLUMNS,transmission_time_gamma psi, transmission_time_gamma* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED,const std::vector<double>& check_times,int seed){

    // Lambda function to convert 
    auto coordinates = [](node_t node, int COLUMNS) -> std::pair<int, int> {
        int quotient = node / COLUMNS;
        int remainder = node % COLUMNS;
        return std::make_pair(quotient, remainder);
    };

    std::vector<double> checking_times;
    for (auto fi : check_times){
        checking_times.push_back(fi);
    }

    int file_number = 0;
    std::string filename = "DATA/snapshot";
    std::string dat = ".dat";
    rng_t engine;


    engine.seed(seed);
    networkx network(graph);
    std::vector<std::vector<int>> grid(ROWS, std::vector<int>(COLUMNS, 0));

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

    std::vector<double> time_trajectory({});
    std::vector<int> infected_trajectory({});

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    // make it random: uncomment
    // for (node_t i = 0; i < INITIAL_INFECTED; i++)
    // {
    //     const node_t random_node = uniform_node_distribution(engine);
    //     // sample without replacement:
    //     if (simulation.is_infected(random_node)){
    //         i--;
    //         continue;
    //     }
    //     simulation.add_infections({ std::make_pair(random_node, 0)});
    //     std::pair<int, int> node_coordinate = coordinates(random_node, COLUMNS);
    //     grid[node_coordinate.first][node_coordinate.second] = 1;
    // }

    const node_t initial_node = std::round((COLUMNS -1)/2 * (COLUMNS + 1));
    std::pair<int, int> node_coordinate = coordinates(initial_node, COLUMNS);
    simulation.add_infections({ std::make_pair(initial_node, 0)});
    grid[node_coordinate.first][node_coordinate.second] = 1;
    
    int current_number_infected = INITIAL_INFECTED - 1;
                        
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;

        const std::pair<int, int> node_coordinate = coordinates(point->node, COLUMNS);
        
        switch (point->kind) {
            case event_kind::outside_infection:
            case event_kind::infection:
                current_number_infected++;
                grid[node_coordinate.first][node_coordinate.second] = 1;
                break;
            case event_kind::reset:
                current_number_infected--;
                grid[node_coordinate.first][node_coordinate.second] = 0;
                break;
            default:
                break;
        }
        if (current_number_infected<=0 || point->time > TMAX){
            py::print("no more infected !");
            break;
        }
        time_trajectory.push_back(point->time);
        infected_trajectory.push_back(current_number_infected);

        // if we are at a checking time, take a snapshot of the epidemic and save it in a dat file.
        double t = std::round(point->time);
        auto it = std::find(checking_times.begin(), checking_times.end(), t);
        if (it != checking_times.end()) {
            checking_times.erase(it);
            file_number ++;
            py::print(file_number,"\r",py::arg("end") = "");
            const std::string temp_filename = filename + std::to_string(file_number) + dat;
            save_grid(grid,temp_filename);
        }



    }
    return std::make_tuple(time_trajectory, infected_trajectory);
}

