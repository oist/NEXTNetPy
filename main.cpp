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

class networkx : public virtual graph_adjacencylist {
public:
    //Initialise graph with adjacency list
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

    int SIZE = py::int_( graph.attr("size")() );
    SIZE +=1;
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
    "A function that runs a SI epidemic on a graph G.");

    py::class_<transmission_time_gamma>(handle, "time_distribution")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0)
        .def_readonly("mean", &transmission_time_gamma::mean)
        .def_readonly("variance", &transmission_time_gamma::variance);


}