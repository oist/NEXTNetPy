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

