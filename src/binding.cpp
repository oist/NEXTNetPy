#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"
#include "dynamic_graph.h"

#include "networkx.hpp"
#include "simulation_wrapper.hpp"
#include "tools.hpp"

namespace py = pybind11;

PYBIND11_MODULE(nmepinet, handle) {

    handle.doc() = "nmepinet module to efficiently simulate an epidemic on any networkx graph, custom graph, or temporal graph.";

    handle.def("simulation_discrete",&simulation_discrete,
        py::arg("graph"),
        py::arg("nb_simulations")=1,
        py::arg("seed")=1,
        py::arg("verbose")=false,
        "returns a vector of vectors (N x kmax) of the number of newly infected with degree k at each discrete step"
    );
    
    handle.def("simulate", 
        py::overload_cast<graph&, transmission_time&, transmission_time*, bool, double, bool, int, int>(&simulate),
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0,
        ""
    );

    handle.def("simulate", 
        py::overload_cast<py::object, transmission_time&, transmission_time*, bool, double, bool, int, int>(&simulate),
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0,
        ""
    );

    handle.def("simulate_average",
        py::overload_cast<graph&, transmission_time&, transmission_time*,bool, double, bool, int, int,int,bool,bool,bool>(&run_simulation_average), 
        py::arg("graph"),
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
        "Simulate average trajectory for custom NMepinet graph object"
    );

    handle.def("simulate_average",
        py::overload_cast<py::object, transmission_time&, transmission_time*, bool, double, bool, int, int,int,bool,bool,bool>(&run_simulation_average), 
        py::arg("graph"),
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
        "Simulate average trajectory on a networkx graph"
    );

    handle.def("simulate_on_temporal", &simulate_on_temporal,
        // py::overload_cast<dynamic_empirical_network&, transmission_time&, transmission_time*, bool, double, bool, int>(&simulate_on_temporal),
        py::arg("temporal_graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("seed")=0,
        py::arg("outside_infection_probability")=0.0,
        ""
    );


    /*---------------------------*/
    /* Graphs that are clustered */
    /*---------------------------*/

    handle.def("graph_ER_clustered",&graph_ER_clustered,
        py::arg("size"),
        py::arg("p"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        ""
    );

    handle.def("graph_LOG_clustered",&graph_LOG_clustered,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        ""
    );

    handle.def("graph_LOG_CM",&graph_LOG_CM,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("r")=0,
        py::arg("seed")=1,
        ""
    );

    handle.def("graph_ER_correlated",&graph_ER_correlated,
        py::arg("size"),
        py::arg("mean"),
        py::arg("r")=0,
        py::arg("seed")=1,
        ""
    );


    //---------------------------------
    //-----------TOOLS ----------------
    //---------------------------------
    
    
    handle.def("degree_clustering_coefficient",&degree_clustering_coefficient,
        py::arg("graph"),
        "computes the average clustering coefficient for nodes of degree k"
    );

    handle.def("connectivity_matrix",&connectivity_matrix,
        py::arg("graph"),
        py::arg("clustering")=0,
        ""
    );

    //---------------------------------
    //---------------------------------
    //---------------------------------


    py::class_<rng_t>(handle, "rng")
        .def(py::init<>()) 
        .def(py::init<int>(),py::arg("seed"),
        ""
        );
      
    py::class_<graph>(handle, "graph", py::multiple_inheritance())
        .def("nodes", &graph::nodes)
        .def("neighbour", &graph::neighbour)
        .def("outdegree", &graph::outdegree);
    
    // py::class_<graph_adjacencylist, graph>(handle, "graph_adjacencylist", py::multiple_inheritance())
    //     .def_readonly("adjacencylist", &graph_adjacencylist::adjacencylist, py::return_value_policy::reference_internal)
    //     .def("al_len", &al_len);    


    // int al_len(graph_adjacencylist* al) {
    //     return al->adjacencylist.size();
    // }


    py::class_<graph_adjacencylist, graph>(handle, "graph_adjacencylist", py::multiple_inheritance())
        .def_readonly("adjacencylist", &graph_adjacencylist::adjacencylist, py::return_value_policy::reference_internal)
        .def("al_len", [](const graph_adjacencylist* al) {
            return al->adjacencylist.size();
        });

    py::class_<watts_strogatz,graph_adjacencylist>(handle, "watts_strogatz", py::multiple_inheritance())
        .def(py::init<node_t, int,double,rng_t&>(), py::arg("size"), py::arg("k"), py::arg("p"),py::arg("rng"),
        "");

    py::class_<erdos_reyni,graph_adjacencylist>(handle, "erdos_reyni", py::multiple_inheritance())
        .def(py::init<int,double,rng_t&>(), py::arg("size"), py::arg("average_degree"),py::arg("rng"),
        "");

    py::class_<barabasi_albert,graph_adjacencylist>(handle, "barabasi_albert", py::multiple_inheritance())
        .def(py::init<int,rng_t&,int>(), py::arg("size"),py::arg("rng"),py::arg("m")=1,
        "");

    py::class_<dynamic_empirical_network, graph>(handle, "temporal_empirical_graph", py::multiple_inheritance())
        .def(py::init<std::string, double>());
        // .def("present_edges", &dynamic_empirical_network::present_edges,
        //     py::arg("t"),
        //     "return number of edges existing at time t"
        // )
        // .def("average_degree", &dynamic_empirical_network::average_degree, 
        //     py::arg("t"), 
        //     "Return average degree of nodes that are active at time t"
        // );

    py::class_<networkx, graph>(handle, "graph_networkx", py::multiple_inheritance())
        .def(py::init<py::list>())
        .def(py::init<py::object>());
    
    py::class_<transmission_time>(handle, "transmission_time")
        .def("sample", &transmission_time::sample);

    py::class_<transmission_time_gamma, transmission_time>(handle, "transmission_time_gamma")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0,
        
        "time_transmission class used to describe the times of infection and/or times of recovery of an individual. \n"
        "For now, the time_distribution is a Gamma distribution by default and only choice."
        "\n"
        "Args:\n"
        "   mean of the distribution.\n"
        "   variance of the distribution.\n"
        "   pinf (double) = 0: probability that the event never happens. For pinf=0 the distribution is well-normalised.\n"
        "However, some functions are not taking this parameter into consideration, for now it is not advised to change pinf.\n"
        )
        .def_readonly("mean", &transmission_time_gamma::mean)
        .def_readonly("variance", &transmission_time_gamma::variance);

    py::class_<transmission_time_deterministic, transmission_time>(handle, "transmission_time_deterministic")
        .def(py::init<double>(), py::arg("tau"),
        // .def(py::init<double,double>(), py::arg("tau"),py::arg("pinf") = 0.0,
        "\n"
        "Args:\n"
        "   deterministic time tau.\n"
        "\n"
        )
        .def_readonly("tau", &transmission_time_deterministic::value);
        // .def_readonly("pinf", &transmission_time_deterministic::pinfinity);


    py::class_<transmission_time_weibull, transmission_time>(handle, "transmission_time_weibull")
        .def(py::init<double, double, double>(), py::arg("shape"), py::arg("scale"), py::arg("pinf") = 0.0,
        
        "time_transmission class used to describe the times of infection and/or times of recovery of an individual. \n"
        "For now, the time_distribution is a Gamma distribution by default and only choice."
        "\n"
        "Args:\n"
        "   mean of the distribution.\n"
        "   variance of the distribution.\n"
        "   pinf (double) = 0: probability that the event never happens. For pinf=0 the distribution is well-normalised.\n"
        "However, some functions are not taking this parameter into consideration, for now it is not advised to change pinf.\n"
        )
        .def_readonly("mean", &transmission_time_weibull::mean)
        .def_readonly("variance", &transmission_time_weibull::variance);

    py::class_<transmission_time_lognormal, transmission_time>(handle, "transmission_time_lognormal")
        .def(py::init<double, double, double>(), py::arg("mean"), py::arg("variance"), py::arg("pinf") = 0.0,
        
        "time_transmission class used to describe the times of infection and/or times of recovery of an individual. \n"
        "For now, the time_distribution is a Gamma distribution by default and only choice."
        "\n"
        "Args:\n"
        "   mean of the distribution.\n"
        "   variance of the distribution.\n"
        "   pinf (double) = 0: probability that the event never happens. For pinf=0 the distribution is well-normalised.\n"
        "However, some functions are not taking this parameter into consideration, for now it is not advised to change pinf.\n"
        )
        .def_readonly("mean", &transmission_time_lognormal::mean)
        .def_readonly("variance", &transmission_time_lognormal::variance);

}
