#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "nextnet/network.h"

namespace py = pybind11;


// Wrapper to convert a python networkx network to our c++ network object.
class networkx : public virtual adjacencylist_network {
public:
    //Initialise network with an adjacency list
    networkx(py::list py_list);

    // OR Initialise network with networkx object
    networkx(py::object network);
};
