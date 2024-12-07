#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "network.h"

namespace py = pybind11;


// Wrapper to convert a python networkx network to our c++ network object.
class networkx : public virtual adjacencylist_network {
public:
    //Initialise network with an adjacency list
    networkx(py::list py_list);

    // OR Initialise network with networkx object
    networkx(py::object network);
};

// Wrapper function to convert a c++ network to a python networkx object
py::object network_ER_clustered(int size,double p,double alpha,double beta,int seed=1);
py::object network_LOG_clustered(int size,double mean,double variance,double alpha,double beta,int seed=1);
py::object network_LOG_CM(int size,double mean,double variance,double r,int seed=1);
py::object network_ER_correlated(int size,double mean,double r,int seed=1);