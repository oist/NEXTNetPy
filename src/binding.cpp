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

    handle.def("simulation_discrete",&simulation_discrete,
        py::arg("network"),
        py::arg("nb_simulations")=1,
        py::arg("seed")=1,
        py::arg("verbose")=false,
        "returns a vector of vectors (N x kmax) of the number of newly infected with degree k at each discrete step"
    );
    
    handle.def("simulate", 
        py::overload_cast<network&, transmission_time&, transmission_time*, bool, double, bool, int, int>(&simulate),
        py::arg("network"),
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
        py::arg("network"),
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
        py::overload_cast<network&, transmission_time&, transmission_time*,bool, double, bool, int, int,int,bool,bool,bool>(&run_simulation_average), 
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
        "Simulate average trajectory for custom nextnet network object"
    );

    handle.def("simulate_average",
        py::overload_cast<py::object, transmission_time&, transmission_time*, bool, double, bool, int, int,int,bool,bool,bool>(&simulate_average), 
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
        "Simulate average trajectory on a networkx network"
    );

    handle.def("simulate_on_activity_driven",&simulate_on_activity_average,
        py::arg("activity_rates"),
        py::arg("inactivation_rate"),
        py::arg("eta"),
        py::arg("m"),
        py::arg("beta"),
        py::arg("mu"),
        py::arg("TMAX")=1000,
        py::arg("seed")=0,
        py::arg("initial_infected")=1,
        py::arg("t0")=0,
        py::arg("nb_simulations")=1,
        py::arg("SIR")=true,
        ""
    );
    // std::tuple<std::vector<double>, std::vector<double>> simulate_on_activity_average(std::vector<double>& activity_rates, double inactivation_rate,double eta, int m, double beta, double mu,double TMAX,int seed,int initial_infected,double t0,int nb_simulations);


    handle.def("simulate_on_temporal",
        py::overload_cast<temporal_network&, transmission_time&, transmission_time*, bool, double, bool, int,int,int,bool,bool,double>(&simulate_on_temporal),
        py::arg("temporal_network"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("seed")=0,
        py::arg("initial_infected")=1,
        py::arg("network_size")=1,
        py::arg("trim")=true,
        py::arg("verbose")=false,
        py::arg("t0")=0,
        ""
    );

    /*---------------------------*/
    /* networks that are clustered */
    /*---------------------------*/

    handle.def("network_ER_clustered",&network_ER_clustered,
        py::arg("size"),
        py::arg("p"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        ""
    );

    handle.def("network_LOG_clustered",&network_LOG_clustered,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        ""
    );

    handle.def("network_LOG_CM",&network_LOG_CM,
        py::arg("size"),
        py::arg("mean"),
        py::arg("variance"),
        py::arg("r")=0,
        py::arg("seed")=1,
        ""
    );

    handle.def("network_ER_correlated",&network_ER_correlated,
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
        py::arg("network"),
        "computes the average clustering coefficient for nodes of degree k"
    );

    handle.def("connectivity_matrix",&connectivity_matrix,
        py::arg("network"),
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
      
    py::class_<network>(handle, "network", py::multiple_inheritance())
        .def("nodes", &network::nodes)
        .def("neighbour", &network::neighbour)
        .def("outdegree", &network::outdegree);

    py::class_<temporal_network>(handle, "temporal_network", py::multiple_inheritance())
        .def("next", &temporal_network::next)
        .def("step", &temporal_network::step)
        .def("notify_epidemic_event", &temporal_network::notify_epidemic_event);

//     struct temporal_network : public virtual network {
// 	virtual absolutetime_t next(rng_t& engine) = 0;

// 	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN) = 0;

// 	virtual void notify_epidemic_event(event_t ev, rng_t& engine);
// };

    
    // py::class_<adjacencylist_network, network>(handle, "adjacencylist_network", py::multiple_inheritance())
    //     .def_readonly("adjacencylist", &adjacencylist_network::adjacencylist, py::return_value_policy::reference_internal)
    //     .def("al_len", &al_len);    


    // int al_len(adjacencylist_network* al) {
    //     return al->adjacencylist.size();
    // }


    py::class_<adjacencylist_network, network>(handle, "adjacencylist_network", py::multiple_inheritance())
        .def_readonly("adjacencylist", &adjacencylist_network::adjacencylist, py::return_value_policy::reference_internal)
        .def("al_len", [](const adjacencylist_network* al) {
            return al->adjacencylist.size();
        });

    py::class_<watts_strogatz,adjacencylist_network>(handle, "watts_strogatz", py::multiple_inheritance())
        .def(py::init<node_t, int,double,rng_t&>(), py::arg("size"), py::arg("k"), py::arg("p"),py::arg("rng"),
        "");

    py::class_<erdos_reyni,adjacencylist_network>(handle, "erdos_reyni", py::multiple_inheritance())
        .def(py::init<int,double,rng_t&>(), py::arg("size"), py::arg("average_degree"),py::arg("rng"),
        "");

    py::class_<barabasi_albert,adjacencylist_network>(handle, "barabasi_albert", py::multiple_inheritance())
        .def(py::init<int,rng_t&,int>(), py::arg("size"),py::arg("rng"),py::arg("m")=1,
        "");


    py::class_<activity_driven_network,temporal_network>(handle, "activity_driven_network", py::multiple_inheritance())
        .def(py::init<std::vector<double>,double,double,double,rng_t&>(), py::arg("activity_rates"),py::arg("eta"),py::arg("m")=1,py::arg("recovery_rate"),py::arg("rng"),
        "initialise instance of an activity-driven network")
        .def("advance_time", 
         &activity_driven_network::advance_time, 
         py::arg("engine"), py::arg("max_time") = std::numeric_limits<double>::quiet_NaN(),
         "Simulate the network until equilibrium in the average degree is reached or the specified max_time.");

    py::class_<temporal_sirx_network,temporal_network>(handle, "temporal_sirx_network", py::multiple_inheritance())
        .def(py::init<network&, double,double>(), py::arg("network"),py::arg("kappa0"),py::arg("kappa"),
        "");

//  temporal_sirx_network : public virtual temporal_network

    py::class_<empirical_temporal_network, temporal_network>(handle, "empirical_temporal_network", py::multiple_inheritance())
        .def(py::init([](std::string path_to_file, bool is_finite_duration, double dt) {
            // Map the bool to the appropriate enum value
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

    py::class_<transmission_time_exponential, transmission_time>(handle, "transmission_time_exponential")
        .def(py::init<double>(), py::arg("rate"),
        ""
        );
}
