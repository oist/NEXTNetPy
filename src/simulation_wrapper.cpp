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

// rng_t engine;
// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>,std::vector<double>> simulation_discrete(py::object graph,int sim, int seed){
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
        
        //initial condition recursion
        // average_leaf_degree[0] += (double) network.outdegree(random_node) / sim ;


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


        // if the epidemic never went exponential (died after a few steps) then skip that step
        // note: this part of the code is never needed in BA networks as the network is connected.
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
    return std::make_tuple(zn_average,recursion);
}
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> depletion_discrete(py::object graph,int sim, int seed){
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
    std::vector<double> k1_traj(n_min,0);
    std::vector<double> k2_traj(n_min,0);
    std::vector<double> k3_traj(n_min,0);
    std::vector<double> r_traj(n_min,0);

    std::vector<double> k1_leaves(n_min,0);
    std::vector<double> k2_leaves(n_min,0);
    std::vector<double> k3_leaves(n_min,0);
    std::vector<double> r_leaves(n_min,0);

    double k2 = 0;
    double k1 = 0;
    double k3 = 0;
    int kmax = 0;
    for (node_t node = 0; node < SIZE; node++ ){
        const int k = network.outdegree(node);
        kmax = std::max(k,kmax);
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3);
    }
    const double original_k1 = k1;
    const double original_k2 = k2;
    const double original_k3 = k3;

    for (int s = 0; s<sim;s++){
        py::print(s,"/",sim,"\r",py::arg("end") = "");
        
        // scale_free nw(SIZE,engine);
        simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
        std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
        const node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});

        // reinitialise moments for each simulation.
        k1 = original_k1 ;
        k2 = original_k2 ;
        k3 = original_k3;
        
        int current_step = 0;
        int current_infected = 0;
        int n = 0;
        int SIZE_LEFT = SIZE;
        // begin simulation
        std::vector<node_t> leaves;
        // const double r = assortativity_depleted(network,&simulation);
        // py::print("debut: ",r,"\n");
        while (true) {
            auto point = simulation.step(engine);
            if (!point)
                break;

            // update current moments
            SIZE_LEFT --;
            const node_t infected_node = point->node;
            const int k = network.outdegree(infected_node);
            k1 -= k ;
            k2 -= pow(k,2) ;
            k3 -= pow(k,3);
            n = (int) std::round(point->time);
            
            // py::print("step: ",n,"  curr: ",current_step,"\n");
            // count if the infected is part of generation n or if we have reached the next gen
            if (n==current_step){
                current_infected ++;
                // add node to the current set of leaves
                leaves.push_back(point->node);
            } else if (n > current_step){
                zn_average[current_step] += (double) current_infected / sim;

                // compute moments
                k1_traj[current_step] += k1 / (sim * SIZE_LEFT);
                k2_traj[current_step] += k2 / (sim * SIZE_LEFT);
                k3_traj[current_step] += k3 / (sim * SIZE_LEFT);

                std::vector<double> result = analyse_leaves(leaves,network,&simulation,kmax);

                // py::print("average degree of leaf : ",result[1],"\n");

                r_leaves[current_step] += result[0] / sim;
                k1_leaves[current_step] += result[1] / sim ;
                k2_leaves[current_step] += result[2] / sim ;
                k3_leaves[current_step] += result[3] / sim ;

                // clear the set of leaves
                leaves.clear();
                leaves.push_back(point->node);
                // py::print("degree of new leaf:", network.outdegree(point->node),"\n");
                // update step
                current_infected = 1;
                current_step = n;
            } else{
                throw std::logic_error("new step cannot be smaller than old step");
            }
        }
        // epidemic ended, but we need to update one last time the trajectories:
        zn_average[current_step] += (double) current_infected / sim;
        k1_traj[current_step] += k1 / (sim * SIZE_LEFT);
        k2_traj[current_step] += k2 / (sim * SIZE_LEFT);
        k3_traj[current_step] += k3 / (sim * SIZE_LEFT);


        std::vector<double> result = analyse_leaves(leaves,network,&simulation,kmax);

        r_leaves[current_step] += result[0] / sim;
        k1_leaves[current_step] += result[1] / sim ;
        k2_leaves[current_step] += result[2] / sim ;
        k3_leaves[current_step] += result[3] / sim ;

        // if the epidemic never went exponential (died after a few steps) then skip that step
        // note: this part of the code is never needed in BA networks as the network is connected.
        if (n<=5)
            continue;
        n_min = std::min(n,n_min);
        // py::print(n_min);
    }
    while (zn_average.size() > n_min ){
        zn_average.pop_back();
        k1_traj.pop_back();
        k2_traj.pop_back();
        k3_traj.pop_back();
        r_traj.pop_back();
        k1_leaves.pop_back();
        k2_leaves.pop_back();
        k3_leaves.pop_back();
        r_leaves.pop_back();
    }
    return std::make_tuple(zn_average,k2_traj,k2_leaves,k1_traj,k1_leaves);
}




// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> run_simulation(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed){

    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

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


    for (node_t i = 0; i < INITIAL_INFECTED; i++)
    {
        const node_t random_node = uniform_node_distribution(engine);
        // sample without replacement:
        if (simulation.is_infected(random_node)){
            i--;
            continue;
        }
        simulation.add_infections({ std::make_pair(random_node, 0)});
    }
    
    int current_number_infected = INITIAL_INFECTED - 1;
                        
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;

        switch (point->kind) {
            case event_kind::outside_infection:
            case event_kind::infection:
                current_number_infected++;
                break;
            case event_kind::reset:
                current_number_infected--;
                break;
            default:
                break;
        }
        if (current_number_infected<=0 || point->time > TMAX)
            break;
        time_trajectory.push_back(point->time);
        infected_trajectory.push_back(current_number_infected);
    }
    return std::make_tuple(time_trajectory, infected_trajectory);
}



// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool TRIM){
    
    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    // Make a vector of pair, which is more convienient once we need to sort the times.
    std::vector<std::pair<double,double>> time_trajectory;
    std::vector<double> infected_trajectory;

    for (int sim = 0; sim < NB_SIMULATIONS; sim++)
        {
            // std::cout << sim << "/" << NB_SIMULATIONS << "\r";
            py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");
            // Initialise/reset simulation object
            simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
            
            for (node_t i = 0; i < INITIAL_INFECTED; i++)
            {
                const node_t random_node = uniform_node_distribution(engine);
                // sample without replacement:
                if (simulation.is_infected(random_node)){
                    i--;
                    continue;
                }
                simulation.add_infections({ std::make_pair(random_node, 0)});
            }
            
            double infection_counter = INITIAL_INFECTED / NB_SIMULATIONS;
                                
            for (int step = 0; step < SIZE*NB_SIMULATIONS;step++) {
                auto point = simulation.step(engine);
                
                if (!point )
                    break;

                switch (point->kind) {
                    case event_kind::outside_infection:
                    case event_kind::infection:
                        infection_counter =  1.0 ;
                        break;
                    case event_kind::reset:
                        infection_counter =  - 1.0 ;
                        break;
                    default:
                        break;
                }
                if (point->time > TMAX)
                    break;
                time_trajectory.push_back(std::make_pair(point->time , infection_counter ) );
                // infected_trajectory.push_back(current_number_infected/NB_SIMULATIONS);
            }
        }


    std::vector<double> time_trajectory_trimmed;
    std::vector<double> infected_trajectory_trimmed;

    const int COUNTER = (NB_SIMULATIONS > 1 && TRIM) ? NB_SIMULATIONS : 1;
    const int NORM = (NB_SIMULATIONS > 1 && TRIM) ? 1 : NB_SIMULATIONS;

    sort(time_trajectory.begin(), time_trajectory.end(),[](const auto& t1, const auto& t2) {
        return t1.first < t2.first;
    });

    for (int index = 0; index < time_trajectory.size(); index += COUNTER){
            time_trajectory_trimmed.push_back(time_trajectory[index].first);
            infected_trajectory_trimmed.push_back(time_trajectory[index].second / NORM);
    }

    // compute the cumulative sum
    std::partial_sum(infected_trajectory_trimmed.begin(), infected_trajectory_trimmed.end(), infected_trajectory_trimmed.begin());

    return std::make_tuple(time_trajectory_trimmed, infected_trajectory_trimmed);
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

