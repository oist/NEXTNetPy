#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "networkx.hpp"
#include "nextnet/weighted_network.h"

namespace py = pybind11;

// Wrapper to convert a python networkx network to our c++ network object.
networkx::networkx(py::list py_list){

    for (auto row : py_list) {
        std::vector<int> cpp_row;
        for (auto element : row) {
            cpp_row.push_back(py::cast<node_t>(element));
        }
        adjacencylist.push_back(cpp_row);
    }

}

    // OR Initialise network with networkx network
networkx::networkx(py::object network){

    int SIZE = py::int_( network.attr("number_of_nodes")() );

    adjacencylist.resize(SIZE);
    for (int i = 0; i < SIZE; i++) {
        py::iterator it = network.attr("neighbors")(i);
        while(it != py::iterator::sentinel()){
            const int j = py::cast<int>(*it);
            adjacencylist[i].push_back(j);
            ++it;
        }
    }
}

// weighted_networkx::weighted_networkx(py::object network) {

//     int SIZE = py::int_(network.attr("number_of_nodes")());

//     adjacencylist.resize(SIZE);
//     for (int i = 0; i < SIZE; ++i) {
//         py::iterator it = network.attr("neighbors")(i);
//         while (it != py::iterator::sentinel()) {
//             const int j = py::cast<int>(*it);

//             // Get edge data dict for (i, j)
//             py::object edge_data = network.attr("get_edge_data")(i, j);

//             // Default weight if not present or no edge data
//             double w = 1.0;

//             if (!edge_data.is_none()) {
//                 // edge_data is a dict; use .get("weight", default)
//                 py::object w_obj =
//                     edge_data.attr("get")("weight", py::float_(1.0));
//                 w = w_obj.cast<double>();
//             }

//             adjacencylist[i].emplace_back(j, w);
//             ++it;
//         }
//     }
// }

weighted_networkx::weighted_networkx(py::object network) {

    const int SIZE = py::int_(network.attr("number_of_nodes")());

    adjacencylist.resize(SIZE);

    // Cache Python callables/keys
    py::object get_edge_data = network.attr("get_edge_data");
    py::object neighbors     = network.attr("neighbors");
    py::str weight_key("weight");

    for (int i = 0; i < SIZE; ++i) {
        py::iterator it = neighbors(i);
        while (it != py::iterator::sentinel()) {
            const int j = py::cast<int>(*it);

            // get_edge_data(i, j) -> dict or None
            py::object edge_data_obj = get_edge_data(i, j);

            double w = 1.0;  // default

            if (!edge_data_obj.is_none()) {
                py::dict edge_data(edge_data_obj);
                if (edge_data.contains(weight_key)) {
                    w = py::cast<double>(edge_data[weight_key]);
                }
            }

            adjacencylist[i].emplace_back(j, w);
            ++it;
        }
    }
}
    

#if 0

py::object network_ER_correlated(int size,double mean,double r,int seed){
    rng_t engine;
    engine.seed(seed);
    // py::print("creating network...\r",py::arg("end") = "");
    erdos_reyni network(size,mean,engine);


    if (r != 0){
        add_correlation(r,network,engine);
    }
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

    // convert to networkx network object
    py::object networkx = py::module::import("networkx");
    py::object G = networkx.attr("network")();
    // Add all edges to the network at once
    G.attr("add_edges_from")(py_edges);

    return G;

}

py::object network_ER_clustered(int size,double p,double alpha,double beta,int seed){
    rng_t engine(seed);

    // Generate a Poisson network with the configuration model
    std::poisson_distribution<> poisson((double) size * p);
    std::vector<int> degrees(size,0);
    std::size_t total_degree = 0;
    for (int i = 0; i < size; i++)
    {
        const int k = poisson(engine);
        degrees[i] = k;
        total_degree += k;
    }

    // make sure the total degree is even, otherwise no network can exist
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

    // convert to networkx network object
    py::object networkx = py::module::import("networkx");
    py::object G = networkx.attr("network")();
    // Add all edges to the network at once
    G.attr("add_edges_from")(py_edges);

    return G;
}

py::object network_LOG_CM(int size,double mean,double variance,double r, int seed){
    rng_t engine;
    engine.seed(seed);
    py::print("creating network...\r",py::arg("end") = "");
    std::vector<int> degrees = lognormal_degree_list(mean,variance,size,engine);
    py::print("CONFIG...\r",py::arg("end") = "");
    config_model network(degrees,engine);

    if (r != 0){
        add_correlation(r,network,engine);
    }
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

    // convert to networkx network object
    py::object networkx = py::module::import("networkx");
    py::object G = networkx.attr("network")();
    // Add all edges to the network at once
    G.attr("add_edges_from")(py_edges);

    return G;

}

py::object network_LOG_clustered(int size,double mean,double variance,double alpha,double beta,int seed){
    rng_t engine(seed);

    std::vector<int> degrees = lognormal_degree_list(mean,variance,size,engine);

 
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

    // convert to networkx network object
    py::object networkx = py::module::import("networkx");
    py::object G = networkx.attr("network")();
    // Add all edges to the network at once
    G.attr("add_edges_from")(py_edges);

    return G;
}

#endif
