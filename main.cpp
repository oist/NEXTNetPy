#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"


namespace py = pybind11;
rng_t engine;

std::vector<std::vector<int>> f(py::list py_list) {
    std::vector<std::vector<int>> cpp_vector;

    for (auto row : py_list) {
        std::vector<int> cpp_row;
        for (auto element : row) {
            cpp_row.push_back(py::cast<int>(element));
        }
        cpp_vector.push_back(cpp_row);
    }

    return cpp_vector;
}

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

void test_function(py::object graph){

    networkx nw(graph);

    py::print(nw.adjacencylist.size());

    for (int i = 0; i < (int) nw.adjacencylist.size(); i++){
        std::cout << "\n";
        for(int v : nw.adjacencylist[i])
            std::cout << v << " ";
    }
    // py::object nn = graph.attr("neighbors");
    // py::iterator it = nn(1);
    // while(it != py::iterator::sentinel()) {
    //     py::print(*it);
    //     ++it;
    // }
}


std::vector<double> run_simulation(py::object graph, int seed = 1){
    
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    const int INITIAL_INFECTED = 1;
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool EDGES_CONCURRENT = true;

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,false);

    std::vector<double> trajectory({});

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);


    for (node_t i = 0; i < INITIAL_INFECTED; i++)
    {
        node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
    }
                        
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;
        trajectory.push_back(point->time);
        
    }

    return trajectory;
}




PYBIND11_MODULE(episimpy, handle) {
    handle.doc() = "episimpy module to efficiently simulate an epidemic on any networkx graph."; // optional module docstring

    handle.def("convert",&f);
    handle.def("simulate", &run_simulation, py::arg("networkx graph"), py::arg("seed")=0, "A function that runs a SI epidemic on a graph G.");

    handle.def("test",&test_function);

    py::class_<transmission_time_gamma>(handle, "time_distribution")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0)
        .def_readonly("mean", &transmission_time_gamma::mean)
        .def_readonly("variance", &transmission_time_gamma::variance);


}