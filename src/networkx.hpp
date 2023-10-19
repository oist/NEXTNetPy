#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "graph.h"

namespace py = pybind11;


// Wrapper to convert a python networkx graph to our c++ graph object.
class networkx : public virtual graph_adjacencylist {
public:
    //Initialise graph with an adjacency list
    networkx(py::list py_list);

    // OR Initialise graph with networkx graph
    networkx(py::object graph);
};

// Wrapper function to convert a c++ graph to a python networkx graph
std::vector<std::vector<size_t>> edges_ER_clustered(int size,double p,double alpha,double beta,int seed=1);
std::vector<std::vector<size_t>> edges_LOG_clustered(int size,double mean,double variance,double alpha,double beta,int seed=1);