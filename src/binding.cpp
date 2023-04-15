#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"

namespace py = pybind11;
rng_t engine;

// Wrapper to convert a python networkx graph to our c++ graph object.
class networkx : public virtual graph_adjacencylist {
public:
    //Initialise graph with an adjacency list
    networkx(py::list py_list){

    for (auto row : py_list) {
        std::vector<int> cpp_row;
        for (auto element : row) {
            cpp_row.push_back(py::cast<node_t>(element));
        }
        adjacencylist.push_back(cpp_row);
    }

    }

    // OR Initialise graph with networkx graph
    networkx(py::object graph){

    int SIZE = py::int_( graph.attr("number_of_nodes")() );
    adjacencylist.resize(SIZE);
    for (int i = 0; i < SIZE; i++) {
        py::iterator it = graph.attr("neighbors")(i);
        while(it != py::iterator::sentinel()){
            const int j = py::cast<int>(*it);
            adjacencylist[i].push_back(j);
            ++it;
        }
    }

    }
};

// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> run_simulation(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1){
    
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
        node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
    }
    
    int current_number_infected = INITIAL_INFECTED;
                        
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
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1, int NB_SIMULATIONS = 1, bool TRIM = false){
    
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    std::vector<double> time_trajectory(SIZE*NB_SIMULATIONS,0);
    std::vector<double> infected_trajectory(SIZE*NB_SIMULATIONS,0);

    for (int sim = 0; sim < NB_SIMULATIONS; sim++)
        {

            // Initialise/reset simulation object
            simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
            
            for (node_t i = 0; i < INITIAL_INFECTED; i++)
            {
                node_t random_node = uniform_node_distribution(engine);
                // sample without replacement:
                if (simulation.is_infected(random_node)){
                    i--;
                    continue;
                }
                simulation.add_infections({ std::make_pair(random_node, 0)});
            }
            
            int64_t current_number_infected = INITIAL_INFECTED;
                                
            for (int step = 0; step < SIZE*NB_SIMULATIONS;step++) {
                auto point = simulation.step(engine);
                if (!point )
                    break;

                switch (point->kind) {
                    case event_kind::outside_infection:
                    case event_kind::infection:
                        current_number_infected+= 1;
                        break;
                    case event_kind::reset:
                        current_number_infected-= 1;
                        break;
                    default:
                        break;
                }
                if (current_number_infected<=0 || point->time > TMAX)
                    break;
                time_trajectory[step] = point->time;
                infected_trajectory[step] = current_number_infected/NB_SIMULATIONS;
            }
        }

    while (time_trajectory.back() <= 0.001){
        time_trajectory.pop_back();
        infected_trajectory.pop_back();
    }

    std::vector<double> time_trajectory_trimmed;
    std::vector<double> infected_trajectory_trimmed;
    if (NB_SIMULATIONS > 1 && TRIM){

        sort(time_trajectory.begin(), time_trajectory.end());

        for (int index = 0; index < time_trajectory.size(); index += NB_SIMULATIONS){
                time_trajectory_trimmed.push_back(time_trajectory[index]);
                infected_trajectory_trimmed.push_back(infected_trajectory[index]);
        }
        return std::make_tuple(time_trajectory_trimmed, infected_trajectory_trimmed);
    } else {
        return std::make_tuple(time_trajectory, infected_trajectory);
    }
}




PYBIND11_MODULE(episimpy, handle) {
    handle.doc() = "episimpy module to efficiently simulate an epidemic on any networkx graph."; // optional module docstring

    handle.def("simulate", &run_simulation, 
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        "Compute the trajectory of a SI/SIR/SIS epidemic on a graph G.");

    handle.def("simulate_average", &run_simulation_average, 
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,

        "Compute the AVERAGE trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
        "\n"
        "Args:\n"
        "    g (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (episimpy.time_distribution): a distribution object that describes the infection times.\n"
        "\n"
        "Returns:\n"
        "    time_traj,infected_traj (list, list): the average trajectory of the epidemic on the graph g.");

    py::class_<transmission_time_gamma>(handle, "time_distribution")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0)
        .def_readonly("mean", &transmission_time_gamma::mean)
        .def_readonly("variance", &transmission_time_gamma::variance);


}