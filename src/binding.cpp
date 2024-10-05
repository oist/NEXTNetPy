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
#include "figures_paper.hpp"
#include "measure_performance.hpp"

namespace py = pybind11;

template <typename T>
void bind_simulation_clustered(py::module &handle) {
    handle.def("simulate_on_clustered",&simulate_on_clustered<T>,
        py::arg("size"),
        py::arg("mean_degree"),
        py::arg("variance"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=1, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,
        py::arg("verbose")=false,
        ""
    );
}


template <typename T>
void bind_simulation_lognormal(py::module &handle) {
    handle.def("simulate_on_lognormal",&simulate_on_lognormal<T>,
        py::arg("size"),
        py::arg("mean_degree"),
        py::arg("variance"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=1, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,
        py::arg("verbose")=false,
        py::arg("r")=0.0,
        ""
    );
}

template <typename T>
void bind_simulation_powerlaw(py::module &handle) {
    handle.def("simulate_on_powerlaw",&simulate_on_powerlaw<T>,
        py::arg("size"),
        py::arg("exponent"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=INFINITY,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=1, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,
        py::arg("add_correlation")=false,
        ""
    );
}


template <typename T>
void bind_simulation_regir(py::module &handle) {
    handle.def("simulate_traj_regir", &simulate_traj_regir<T>, 
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        "Compute a trajectory using REGIR of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
        "\n"
        "Args:\n"
        "    g (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (nmepinet.time_distribution): a distribution object that describes the (attempted) infection times.\n"
        "   recovery_time (nmepinet.time_distribution) = None. a dist. obj. that describes the times of recovery. \n"
        "   SIR (bool) = False. If set to true, the model is with no reinfection (SIR), if false, reinfections are allowed and the simulation describes a SIS model. A recovery_time has to be defined in order for the SIR bool to have an effect. \n"
        "   TMAX (float): value that determines the maximum time an epidemic is being simulated. By default is set to 1000. Even for large graphs, this is very long. We advise to reduce this value if one wants to simulate a SIS epidemic."
        "   concurrent_edges (bool) = False. \n"
        "   initial_infected (int) = 1. Initial number of infected when the epidemic starts at t=0. These nodes are picked uniformly at random among all the nodes of the networks.\n"
        "   seed (int)=0. random seed.\n"
        "   nb_simulations (int) = 1. number of simulations. The graph stays the same, but the epidemic is reset and starts with a different set of infected nodes.\n"
        "  trim (bool)=True. To compute the average. ALL the trajectories are joined and sorted. For large number of simulations, this return very large data files (length SIZE*NB_SIMULATION) if set to False. \n"
        "\n"
        "Returns:\n"
        "    time_traj (list),infected_traj (list): the average trajectory of the epidemic on the graph g.\n"
        
        
    );
}


// template <typename T>
// void bind_simulation(py::module &handle) {
//     handle.def("simulate", &simulate<T>, 
//         py::arg("graph"),
//         py::arg("infection_time"),
//         py::arg("recovery_time")=nullptr,
//         py::arg("SIR")=false,
//         py::arg("TMAX")=1000,
//         py::arg("concurrent_edges")=true,
//         py::arg("initial_infected")=1,
//         py::arg("seed")=0, 
//         "Compute a trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
//         "\n"
//         "Args:\n"
//         "    g (nx.graph): graph generated by the Networkx module.\n"
//         "    infection_time (nmepinet.time_distribution): a distribution object that describes the (attempted) infection times.\n"
//         "   recovery_time (nmepinet.time_distribution) = None. a dist. obj. that describes the times of recovery. \n"
//         "   SIR (bool) = False. If set to true, the model is with no reinfection (SIR), if false, reinfections are allowed and the simulation describes a SIS model. A recovery_time has to be defined in order for the SIR bool to have an effect. \n"
//         "   TMAX (float): value that determines the maximum time an epidemic is being simulated. By default is set to 1000. Even for large graphs, this is very long. We advise to reduce this value if one wants to simulate a SIS epidemic."
//         "   concurrent_edges (bool) = False. \n"
//         "   initial_infected (int) = 1. Initial number of infected when the epidemic starts at t=0. These nodes are picked uniformly at random among all the nodes of the networks.\n"
//         "   seed (int)=0. random seed.\n"
//         "   nb_simulations (int) = 1. number of simulations. The graph stays the same, but the epidemic is reset and starts with a different set of infected nodes.\n"
//         "  trim (bool)=True. To compute the average. ALL the trajectories are joined and sorted. For large number of simulations, this return very large data files (length SIZE*NB_SIMULATION) if set to False. \n"
//         "\n"
//         "Returns:\n"
//         "    time_traj (list),infected_traj (list): the average trajectory of the epidemic on the graph g.\n"
        
        
//     );
// }

// template <typename T>
// void bind_simulation_average(py::module &handle) {
//     handle.def("simulate_average", &run_simulation_average<T>, 
//         py::arg("graph"),
//         py::arg("infection_time"),
//         py::arg("recovery_time")=nullptr,
//         py::arg("SIR")=true,
//         py::arg("TMAX")=1000,
//         py::arg("concurrent_edges")=true,
//         py::arg("initial_infected")=1,
//         py::arg("seed")=0, 
//         py::arg("nb_simulations")=1,
//         py::arg("trim")=true,
//         py::arg("verbose")=false,
//         py::arg("all_nodes")=false,
//         ""
//     );
// }

template <typename T>
void bind_simulation_per_degree_class(py::module &handle) {
    handle.def("simulate_per_degree_class", &simulation_per_degree_class<T>, 
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
        ""
    );
}


template <typename T>
void bind_measure_run_time(py::module &handle) {
    handle.def("measure_run_time",&measure_run_time<T>,
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("verbose")=false,
        ""
    );
}

template <typename T>
void bind_measure_run_time_nmga(py::module &handle) {
    handle.def("measure_run_time_nmga",&measure_run_time_nmga<T>,
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("verbose")=false,
        ""
    );
}

template <typename T>
void bind_measure_run_time_regir(py::module &handle) {
    handle.def("measure_run_time_regir",&measure_run_time_regir<T>,
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("verbose")=false,
        py::arg("approximation")=1,
        ""
    );
}


    int al_len(graph_adjacencylist* al) {
        return al->adjacencylist.size();
    }


PYBIND11_MODULE(nmepinet, handle) {

    handle.doc() = "nmepinet module to efficiently simulate an epidemic on any networkx graph."; // optional module docstring

    // bind_simulation<transmission_time_gamma>(handle);
    // bind_simulation<transmission_time_weibull>(handle);
    // bind_simulation<transmission_time_lognormal>(handle);
    // bind_simulation<transmission_time_exponential>(handle);
    // bind_simulation<transmission_time_deterministic>(handle);

    bind_simulation_regir<transmission_time_gamma>(handle);
    bind_simulation_regir<transmission_time_weibull>(handle);
    bind_simulation_regir<transmission_time_lognormal>(handle);
    bind_simulation_regir<transmission_time_exponential>(handle);
    bind_simulation_regir<transmission_time_deterministic>(handle);

    bind_simulation_per_degree_class<transmission_time_gamma>(handle);
    bind_simulation_per_degree_class<transmission_time_weibull>(handle);
    bind_simulation_per_degree_class<transmission_time_lognormal>(handle);
    bind_simulation_per_degree_class<transmission_time_exponential>(handle);
    bind_simulation_per_degree_class<transmission_time_deterministic>(handle);

    // bind_simulation_average<transmission_time_gamma>(handle);
    // bind_simulation_average<transmission_time_weibull>(handle);
    // bind_simulation_average<transmission_time_lognormal>(handle);
    // bind_simulation_average<transmission_time_exponential>(handle);
    // bind_simulation_average<transmission_time_deterministic>(handle);

    bind_simulation_clustered<transmission_time_gamma>(handle);
    bind_simulation_clustered<transmission_time_weibull>(handle);
    bind_simulation_clustered<transmission_time_lognormal>(handle);
    bind_simulation_clustered<transmission_time_exponential>(handle);
    bind_simulation_clustered<transmission_time_deterministic>(handle);

    bind_simulation_lognormal<transmission_time_gamma>(handle);
    bind_simulation_lognormal<transmission_time_weibull>(handle);
    bind_simulation_lognormal<transmission_time_lognormal>(handle);
    bind_simulation_lognormal<transmission_time_exponential>(handle);
    bind_simulation_lognormal<transmission_time_deterministic>(handle);

    bind_simulation_powerlaw<transmission_time_gamma>(handle);
    bind_simulation_powerlaw<transmission_time_weibull>(handle);
    bind_simulation_powerlaw<transmission_time_lognormal>(handle);
    bind_simulation_powerlaw<transmission_time_exponential>(handle);
    bind_simulation_powerlaw<transmission_time_deterministic>(handle);

    bind_measure_run_time<transmission_time_gamma>(handle);
    bind_measure_run_time<transmission_time_weibull>(handle);
    bind_measure_run_time<transmission_time_lognormal>(handle);
    bind_measure_run_time<transmission_time_exponential>(handle);
    bind_measure_run_time<transmission_time_deterministic>(handle);

    bind_measure_run_time_nmga<transmission_time_gamma>(handle);
    bind_measure_run_time_nmga<transmission_time_weibull>(handle);
    bind_measure_run_time_nmga<transmission_time_lognormal>(handle);
    bind_measure_run_time_nmga<transmission_time_exponential>(handle);
    bind_measure_run_time_nmga<transmission_time_deterministic>(handle);

    bind_measure_run_time_regir<transmission_time_gamma>(handle);
    bind_measure_run_time_regir<transmission_time_weibull>(handle);
    bind_measure_run_time_regir<transmission_time_lognormal>(handle);
    bind_measure_run_time_regir<transmission_time_exponential>(handle);
    bind_measure_run_time_regir<transmission_time_deterministic>(handle);

    handle.def("graph_ER_clustered",&graph_ER_clustered,
        py::arg("size"),
        py::arg("p"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("seed")=1,
        ""
    );

    handle.def("degree_clustering_coefficient",&degree_clustering_coefficient,
        py::arg("graph"),
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

    handle.def("connectivity_matrix",&connectivity_matrix,
        py::arg("graph"),
        py::arg("clustering")=0,
        ""
    );


    handle.def("mu",&mu,
        py::arg("graph"),
        py::arg("clustering")=false,
        ""
    );

    handle.def("depleted_distribution",&depleted_distribution,
        py::arg("SIZE"),
        py::arg("SIM"),
        py::arg("seed"),
        ""
    );

    handle.def("simulation_discrete",&simulation_discrete,
        py::arg("graph"),
        py::arg("nb_simulations")=1,
        py::arg("seed")=1,
        py::arg("verbose")=false,
        "returns a vector of vectors (N x kmax) of the number of newly infected with degree k at each discrete step"
    );

    // handle.def("run_benchmark_next_reaction", &run_benchmark_next_reaction,
    //     py::arg("graph_ensemble"),
    //     py::arg("infection_time"),
    //     py::arg("recovery_time")=nullptr,
    //     py::arg("SIR")=false,
    //     py::arg("TMAX")=1000,
    //     py::arg("concurrent_edges")=true,
    //     py::arg("initial_infected")=1,
    //     py::arg("seed")=0, 
    //     " Average run time of the next reaction method \n"
    //     "\n"
    //     "Returns:\n"
    //     "   list of network sizes, list average run time, (list) std deviation run time"
    // );

    handle.def("generalised_knn",&generalised_knn,
        py::arg("size"),
        py::arg("sim")=1,
        py::arg("power")=2,
        py::arg("seed")=1,
        " nth moment of the degree distribution of the NEIGHBOURS of a node of degree k=0,1,...kmax"
        );

    handle.def("m_nn",&neighbours_multiplicity,
        py::arg("graph"),
        ""
        );
    
    handle.def("edge_multiplicity",&edge_multiplicity2,
        py::arg("graph"),
        ""
        );

    handle.def("simulate_discrete_leaves",&simulation_discrete_leaves,
        py::arg("graph"),
        py::arg("nb_sim")=100,
        py::arg("seed")=0,
        ""
    );

    // handle.def(assortativity_depleted),&assortativity_depleted,
    //     p

    handle.def("depletion",&depletion,
        py::arg("graph"),
        py::arg("frequency_list") = std::vector<double>{0.1,0.25,0.5,0.75,0.9},
        py::arg("seed")=0,
        "Function that returns the empirical degree distribution of the susceptible nodes at different stages of the epidemic.\n"
        "\n"
        "Returns:\n"
        "   vector of size (kmax+1,10) that contains prob. deg. distr. of the sus. nodes at various stages of the epidemic.\n"
        "   list of 4 trajectories of the first 4 evolving moments of the sus.nodes.\n"
    "\n");



    handle.def("growth_rate",&euler_lotka_growth_rate,
        py::arg("graph"),
        py::arg("infection_time"),
        "Computes the growth rate that solves the generalised Euler-Lotka equation for a SI epidemic with a gamma distribution.\n"
        "\n"
        "Args:\n"
        "    graph (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (nmepinet.time_distribution): a gamma distribution object that describes the infection times.\n"
        "\n"
        "Returns:\n"
        "    growth_rate (float): the predicted growth rate of the epidemic."
    );

    handle.def("simulate_on_lattice", &run_simulation_lattice, 
        py::arg("graph"),
        py::arg("rows"),
        py::arg("columns"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("checking_times")=std::vector<double>{0,10,20,30,40,50,60,70,80,90,100},
        py::arg("seed")=0, 
        "Compute a trajectory of a SI/SIR/SIS epidemic on a lattice graph G and save snapshots of the epidemic.\n"
        "\n"
        "Args:\n"
        "    g (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (nmepinet.time_distribution): a distribution object that describes the (attempted) infection times.\n"
        "   recovery_time (nmepinet.time_distribution) = None. a dist. obj. that describes the times of recovery. \n"
        "   SIR (bool) = False. If set to true, the model is with no reinfection (SIR), if false, reinfections are allowed and the simulation describes a SIS model. A recovery_time has to be defined in order for the SIR bool to have an effect. \n"
        "   TMAX (float): value that determines the maximum time an epidemic is being simulated. By default is set to 1000. Even for large graphs, this is very long. We advise to reduce this value if one wants to simulate a SIS epidemic."
        "   concurrent_edges (bool) = False. \n"
        "   initial_infected (int) = 1. Initial number of infected when the epidemic starts at t=0. These nodes are picked uniformly at random among all the nodes of the networks.\n"
        "   seed (int)=0. random seed.\n"
        "   nb_simulations (int) = 1. number of simulations. The graph stays the same, but the epidemic is reset and starts with a different set of infected nodes.\n"
        "  trim (bool)=True. To compute the average. ALL the trajectories are joined and sorted. For large number of simulations, this return very large data files (length SIZE*NB_SIMULATION) if set to False. \n"
        "\n"
        "Returns:\n"
        "    time_traj (list),infected_traj (list): the average trajectory of the epidemic on the graph g.\n"
        
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
        ""
    );
// py::object graph,transmission_time& psi, transmission_time* rho, bool SIR,double TMAX, bool EDGES_CONCURRENT,int INITIAL_INFECTED, int seed, int NB_SIMULATIONS, bool 
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
        ""
    );
    // handle.def("simulate", &simulate,
    //     py::arg("graph"),
    //     py::arg("infection_time"),
    //     py::arg("recovery_time")=nullptr,
    //     py::arg("SIR")=false,
    //     py::arg("TMAX")=1000,
    //     py::arg("concurrent_edges")=true,
    //     py::arg("initial_infected")=1,
    //     py::arg("seed")=0,
    //     ""
    // );

    handle.def("simulate_on_temporal", &simulate_on_temporal,
        // py::overload_cast<dynamic_empirical_network&, transmission_time&, transmission_time*, bool, double, bool, int>(&simulate_on_temporal),
        py::arg("temporal_graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=true,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("seed")=0,
        ""
    );

    py::class_<rng_t>(handle, "rng")
        .def(py::init<>()) 
        .def(py::init<int>(),py::arg("seed"),
        ""
        );
      
    py::class_<graph>(handle, "graph", py::multiple_inheritance())
        .def("nodes", &graph::nodes)
        .def("neighbour", &graph::neighbour)
        .def("outdegree", &graph::outdegree);
    
    py::class_<graph_adjacencylist, graph>(handle, "graph_adjacencylist", py::multiple_inheritance())
        .def_readonly("adjacencylist", &graph_adjacencylist::adjacencylist, py::return_value_policy::reference_internal)
        .def("al_len", &al_len);    

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
        .def(py::init<std::string, double>())
        .def("present_edges", &dynamic_empirical_network::present_edges,
            py::arg("t"),
            "return number of edges existing at time t"
        )
        .def("average_degree", &dynamic_empirical_network::average_degree, 
            py::arg("t"), 
            "Return average degree of nodes that are active at time t"
        );

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
