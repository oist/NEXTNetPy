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
py::object graph_ER_clustered(int size,double p,double alpha,double beta,int seed=1);
py::object graph_LOG_clustered(int size,double mean,double variance,double alpha,double beta,int seed=1);
py::object graph_LOG_CM(int size,double mean,double variance,double r,int seed=1);
py::object graph_ER_correlated(int size,double mean,double r,int seed=1);