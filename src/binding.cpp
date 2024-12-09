#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"
#include "temporal_network.h"

#include "networkx.hpp"
#include "simulation_wrapper.hpp"
#include "tools.hpp"

namespace py = pybind11;

PYBIND11_MODULE(nextnet, handle) {

    handle.doc() = "nextnet module to efficiently simulate an epidemic on any networkx network, custom network, or temporal network.";

    handle.def("simulate", 
        py::overload_cast<network&, transmission_time&, transmission_time*, bool, double, bool, int, int>(&simulate),
        py::arg("network"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("tmax")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0,
        R"(
        Simulate the epidemic spread on a given network.

        This function models the spread of an epidemic on a network with the specified parameters.
        It supports various inbuilt network structures. The simulation runs until the maximum time 
        (`tmax`) is reached, or until the epidemic process is completed.
        The function returns two lists:
            1. A list of infection times for each infected individual.
            2. A list of the number of infected individuals at each time step.

        Args:
            network (network): The network structure on which the simulation will run.
            infection_time (transmission_time): The time an individual stays infected before recovering.
            recovery_time (transmission_time*, optional): The time an individual stays in recovery (default: NULL).
            SIR (bool, optional): If TRUE, the SIR model is used; otherwise, a custom model is used (default: TRUE).
            tmax (double, optional): The maximum simulation time (default: 1000).
            concurrent_edges (bool, optional): If TRUE, concurrent edges are allowed in the simulation (default: TRUE).
            initial_infected (int, optional): The number of initially infected individuals (default: 1).
            seed (int, optional): The random seed to ensure reproducibility (default: 0).

        Returns:
            list: A list of infection times for each infected individual.
            list: A list of the number of infected individuals at each time step.
        )"
    );

    handle.def("simulate", 
        py::overload_cast<py::object, transmission_time&, transmission_time*, bool, double, bool, int, int>(&simulate),
        py::arg("network"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0,
        R"(
        Simulate the epidemic spread on a NetworkX-based network.

        This function models the spread of an epidemic on a NetworkX network with the specified parameters.
        The simulation runs until the maximum time (`tmax`) is reached, or until the epidemic process is completed.
        The function returns two lists:
            1. A list of infection times for each infected individual.
            2. A list of the number of infected individuals at each time step.

        Args:
            network (network): The NetworkX structure on which the simulation will run.
            infection_time (transmission_time): The time an individual stays infected before recovering.
            recovery_time (transmission_time*, optional): The time an individual stays in recovery (default: NULL).
            SIR (bool, optional): If TRUE, the SIR model is used; otherwise, a custom model is used (default: TRUE).
            tmax (double, optional): The maximum simulation time (default: 1000).
            concurrent_edges (bool, optional): If TRUE, concurrent edges are allowed in the simulation (default: TRUE).
            initial_infected (int, optional): The number of initially infected individuals (default: 1).
            seed (int, optional): The random seed to ensure reproducibility (default: 0).

        Returns:
            list: A list of infection times for each infected individual.
            list: A list of the number of infected individuals at each time step.
        )"
    );

    handle.def("simulate_average",
        py::overload_cast<py::object, transmission_time&, transmission_time*, bool, double, bool, int, int, int, bool, bool, bool>(&simulate_average), 
        py::arg("network"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,
        py::arg("verbose")=false,
        py::arg("all_nodes")=false,
        R"(
        Simulate average trajectory on a NetworkX-based network.

        This function models the average trajectory of an epidemic spread over multiple simulations. It returns
        the average behavior of infected individuals over several runs.

        Args:
            network (network): The NetworkX structure on which the simulation will run.
            infection_time (transmission_time): The time an individual stays infected before recovering.
            recovery_time (transmission_time*, optional): The time an individual stays in recovery (default: NULL).
            SIR (bool, optional): If TRUE, the SIR model is used; otherwise, a custom model is used (default: TRUE).
            TMAX (double, optional): The maximum simulation time (default: 1000).
            concurrent_edges (bool, optional): If TRUE, concurrent edges are allowed in the simulation (default: TRUE).
            initial_infected (int, optional): The number of initially infected individuals (default: 1).
            seed (int, optional): The random seed to ensure reproducibility (default: 0).
            nb_simulations (int, optional): The number of simulations to run (default: 1).
            trim (bool, optional): If TRUE, trims results to a given range (default: TRUE).
            verbose (bool, optional): If TRUE, provides detailed output (default: FALSE).
            all_nodes (bool, optional): If TRUE, includes all nodes in the result (default: FALSE).

        Returns:
            list: A list of average infection times for each infected individual.
            list: A list of the average number of infected individuals at each time step.
        )"
    );

    handle.def("network_ER_clustered",&network_ER_clustered,
        py::arg("size"),
        py::arg("p"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        R"(
        Generate a clustered version of an Erdős–Rényi network.

        Args:
            size (int): Number of nodes in the network.
            p (double): Probability of edge creation between nodes.
            alpha (double): A parameter controlling the clustering: ck ~ k**-alpha.
            beta (double): A parameter linked with the assortativity of the network.
            seed (int, optional): The random seed to ensure reproducibility (default: 1).

        Returns:
            network: A clustered Erdős–Rényi network.
        )"
    );

    handle.def("network_LOG_clustered",&network_LOG_clustered,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        R"(
        Generate a clustered version of a network whose degrees follow a lognormal distribution.

        Args:
            size (int): Number of nodes in the network.
            mean (double): The mean of the distribution.
            variance (double): The variance of the distribution.
            alpha (double): A parameter controlling the clustering: ck ~ k**-alpha.
            beta (double): A parameter linked with the assortativity of the network.
            seed (int, optional): The random seed to ensure reproducibility (default: 1).

        Returns:
            network: A clustered logarithmic network.
        )"
    );

    handle.def("network_LOG_CM",&network_LOG_CM,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("r")=0,
        py::arg("seed")=1,
        R"(
        Generate a clustered logarithmic network with a CM degree distribution.

        Args:
            size (int): Number of nodes in the network.
            mean (double): The mean of the distribution.
            variance (double): The variance of the distribution.
            r (double, optional): Correlation coefficient (default: 0).
            seed (int, optional): The random seed to ensure reproducibility (default: 1).

        Returns:
            network: A clustered logarithmic network with a CM degree distribution.
        )"
    );

    handle.def("network_ER_correlated",&network_ER_correlated,
        py::arg("size"),
        py::arg("mean"),
        py::arg("r")=0,
        py::arg("seed")=1,
        R"(
        Generate an Erdős–Rényi network with degree correlation r.

        Args:
            size (int): Number of nodes in the network.
            mean (double): The average degree of the nodes.
            r (double, optional): Correlation coefficient (default: 0).
            seed (int, optional): The random seed to ensure reproducibility (default: 1).

        Returns:
            network: A correlated Erdős–Rényi network.
        )"
    );

    //---------------------------------
    //-----------TOOLS ----------------
    //---------------------------------

    // Define degree clustering coefficient function
    handle.def("degree_clustering_coefficient",&degree_clustering_coefficient,
        py::arg("network"),
        R"(
        Computes the average clustering coefficient for nodes of degree k.

        Args:
            network (network): The network structure on which to compute the clustering coefficient.

        Returns:
            double: The average clustering coefficient.
        )"
    );

    // Define connectivity matrix function
    handle.def("connectivity_matrix",&connectivity_matrix,
        py::arg("network"),
        py::arg("clustering")=0,
        R"(
        Computes the connectivity matrix for the network.

        Args:
            network (network): The network structure on which to compute the connectivity matrix.
            clustering (int, optional): Clustering parameter (default: 0).

        Returns:
            matrix: The connectivity matrix of the network.
        )"
    );

    //---------------------------------
    //---------------------------------
    //---------------------------------

    py::class_<rng_t>(handle, "rng")
        .def(py::init<>()) 
        .def(py::init<int>(),py::arg("seed"),
        R"(
        Random number generator class.

        Args:
            seed (int): The random seed to initialize the generator.
        )"
        );
      
    py::class_<network>(handle, "network", py::multiple_inheritance())
        .def("nodes", &network::nodes)
        .def("neighbour", &network::neighbour)
        .def("outdegree", &network::outdegree);

    py::class_<temporal_network>(handle, "temporal_network", py::multiple_inheritance())
        .def("next", &temporal_network::next)
        .def("step", &temporal_network::step)
        .def("notify_epidemic_event", &temporal_network::notify_epidemic_event);

    py::class_<adjacencylist_network, network>(handle, "adjacencylist_network", py::multiple_inheritance())
        .def_readonly("adjacencylist", &adjacencylist_network::adjacencylist, py::return_value_policy::reference_internal)
        .def("al_len", [](const adjacencylist_network* al) {
            return al->adjacencylist.size();
        });

    py::class_<watts_strogatz,adjacencylist_network>(handle, "watts_strogatz", py::multiple_inheritance())
        .def(py::init<node_t, int, double, rng_t&>(), py::arg("size"), py::arg("k"), py::arg("p"), py::arg("rng"));

    py::class_<erdos_reyni,adjacencylist_network>(handle, "erdos_reyni", py::multiple_inheritance())
        .def(py::init<int, double, rng_t&>(), py::arg("size"), py::arg("average_degree"), py::arg("rng"));

    py::class_<barabasi_albert,adjacencylist_network>(handle, "barabasi_albert", py::multiple_inheritance())
        .def(py::init<int, rng_t&, int>(), py::arg("size"), py::arg("rng"), py::arg("m")=1);

    py::class_<activity_driven_network,temporal_network>(handle, "activity_driven_network", py::multiple_inheritance())
        .def(py::init<std::vector<double>, double, double, double, rng_t&>(), py::arg("activity_rates"), py::arg("eta"), py::arg("m")=1, py::arg("recovery_rate"), py::arg("rng"));

    py::class_<temporal_sirx_network,temporal_network>(handle, "temporal_sirx_network", py::multiple_inheritance())
        .def(py::init<network&, double, double>(), py::arg("network"), py::arg("kappa0"), py::arg("kappa"));

    py::class_<empirical_temporal_network, temporal_network>(handle, "empirical_temporal_network", py::multiple_inheritance())
        .def(py::init([](std::string path_to_file, bool is_finite_duration, double dt) {
            empirical_temporal_network::edge_duration_kind contact_type = 
                is_finite_duration ? empirical_temporal_network::finite_duration : empirical_temporal_network::infitesimal_duration;
            return new empirical_temporal_network(path_to_file, contact_type, dt);
        }));

    py::class_<networkx, network>(handle, "network_networkx", py::multiple_inheritance())
        .def(py::init<py::list>())
        .def(py::init<py::object>());
    
    py::class_<transmission_time>(handle, "transmission_time")
        .def("sample", &transmission_time::sample);

    py::class_<transmission_time_gamma, transmission_time>(handle, "transmission_time_gamma")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0);

    py::class_<transmission_time_deterministic, transmission_time>(handle, "transmission_time_deterministic")
        .def(py::init<double>(), py::arg("tau"));

    py::class_<transmission_time_weibull, transmission_time>(handle, "transmission_time_weibull")
        .def(py::init<double, double, double>(), py::arg("shape"), py::arg("scale"), py::arg("pinf") = 0.0);

    py::class_<transmission_time_lognormal, transmission_time>(handle, "transmission_time_lognormal")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0);

    py::class_<transmission_time_exponential, transmission_time>(handle, "transmission_time_exponential")
        .def(py::init<double>(), py::arg("rate"));
}
