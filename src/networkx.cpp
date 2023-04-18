#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "networkx.hpp"

namespace py = pybind11;

// Wrapper to convert a python networkx graph to our c++ graph object.
networkx::networkx(py::list py_list){

    for (auto row : py_list) {
        std::vector<int> cpp_row;
        for (auto element : row) {
            cpp_row.push_back(py::cast<node_t>(element));
        }
        adjacencylist.push_back(cpp_row);
    }

}

    // OR Initialise graph with networkx graph
networkx::networkx(py::object graph){

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
