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


py::object graph_ER_clustered(int size,double p,double alpha,double beta,int seed){
    rng_t engine(seed);

    // Generate a Poisson graph with the configuration model
    std::poisson_distribution<> poisson((double) size * p);
    std::vector<int> degrees(size,0);
    std::size_t total_degree = 0;
    for (int i = 0; i < size; i++)
    {
        const int k = poisson(engine);
        degrees[i] = k;
        total_degree += k;
    }

    // make sure the total degree is even, otherwise no graph can exist
    while (total_degree % 2 == 1) {
        // re-generate a random degree
        const std::size_t i = std::uniform_int_distribution<>(0, size-1)(engine);
        const int d = degrees[i];
        const int dp = poisson(engine);
        degrees[i] = dp;
        total_degree += dp - d;
    }

    config_model_clustered_serrano network(degrees,alpha,beta,engine);

    // create edge list
    py::list py_edges;
    // std::vector<std::vector<size_t>> edges;
    // edges.reserve(total_degree/2);

    const std::size_t n = network.adjacencylist.size();
	for(std::size_t i=0; i < n; ++i) {
		const std::unordered_set<node_t> neighbours(network.adjacencylist[i].begin(), network.adjacencylist[i].end());
		for(std::size_t j=0; j < i; ++j) {
			if (neighbours.find(j) != neighbours.end()){
                const py::tuple py_edge = py::make_tuple(i,j);
                py_edges.append(py_edge);
            }
		}
	}


    // return edges;

    // convert to networkx graph object
    py::object networkx = py::module::import("networkx");
    py::object G = networkx.attr("Graph")();
    // Add all edges to the graph at once
    G.attr("add_edges_from")(py_edges);

    return G;
}


std::vector<std::vector<size_t>> edges_LOG_clustered(int size,double mean,double variance,double alpha,double beta,int seed){
    rng_t engine(seed);

    std::vector<int> degrees = lognormal_degree_list(mean,variance,size,engine);

    config_model_clustered_serrano network(degrees,alpha,beta,engine);

    // create edge list
    std::vector<std::vector<size_t>> edges;
    int total_degree = std::accumulate(degrees.begin(), degrees.end(), 0);
    edges.reserve(total_degree/2);

    const std::size_t n = network.adjacencylist.size();
	for(std::size_t i=0; i < n; ++i) {
		const std::unordered_set<node_t> neighbours(network.adjacencylist[i].begin(), network.adjacencylist[i].end());
		for(std::size_t j=0; j < i; ++j) {
			if (neighbours.find(j) != neighbours.end())
                edges.push_back(std::vector<size_t> {i,j});
		}
	}
    return edges;

    // // convert to networkx graph object
    // py::object networkx = py::module::import("networkx");
    // py::object G = networkx.attr("Graph")();
    // // Add all edges to the graph at once
    // G.attr("add_edges_from")(edges);

    // return G;
}