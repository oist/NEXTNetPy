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

std::tuple<std::vector<std::vector<double>>,std::vector<double>,std::vector<double>> simulation_discrete_zk(py::object graph,int SIZE,int sim, int seed){
    rng_t engine;
    engine.seed(seed);

    // networkx network(graph);;
    // const int mean = 5.5 ;
    // const int variance = 60;
    // std::vector<int> degrees = lognormal_degree_list(mean,variance,SIZE,engine);

    // config_model network(degrees,engine);
    // add_correlation(-0.4,network,engine);
    // scale_free network(SIZE,engine);

    erdos_reyni network(SIZE,6,engine);
    add_correlation(0.4,network,engine);
    
    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double r = assortativity(network);
    
    int kmax = 0;
    std::vector<double> pk(SIZE,0);

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

    
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    transmission_time_deterministic psi(1);
    
    int n_min = 50;
    std::vector<std::vector<double>> zn_average(n_min,std::vector<double>(kmax+1,0));
    std::vector<double> step_counter(n_min,0);
    for (int s = 0; s<sim;s++){
        py::print(s,"/",sim,"\r",py::arg("end") = "");
        
        simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
        std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
       
        const node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
        
        
        int n = 0;
        int step = 0;
        // begin simulation
        while (true) {
            auto point = simulation.step(engine);
            
            if (!point)
                break;
            
            const node_t infected_node = point->node;
            const int k = network.outdegree(infected_node);
            n = (int) std::round(point->time);
            
            zn_average[n][k] ++ ;
           
            if (n>step){
                step_counter[step]++;
                step = n;
            }
        }
        step_counter[n]++;

    }

    //delete unused steps
    while (step_counter.back()==0 ){
        step_counter.pop_back();
        zn_average.pop_back();
    }

    //normalise zn:
    for (int n = 0; n<step_counter.size(); n++){
        const double norm = (double) step_counter[n];
        for( int k =0;k<= kmax;k++){
            zn_average[n][k] /= norm;
        }
    }

    return std::make_tuple(zn_average,std::vector<double>{k1,k2,k3,r},pk);
}




// rng_t engine;
// Wrapper to run a simulation given a network and transmission distribution
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

