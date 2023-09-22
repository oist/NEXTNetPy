#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"
#include "networkx.hpp"
#include "simulation_wrapper.hpp"
#include "tools.hpp"
#include "figures_paper.hpp"

namespace py = pybind11;

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
        ""
    );
}

template <typename T>
void bind_simulation(py::module &handle) {
    handle.def("simulate", &simulate<T>, 
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        "Compute a trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
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

PYBIND11_MODULE(nmepinet, handle) {

    handle.doc() = "nmepinet module to efficiently simulate an epidemic on any networkx graph."; // optional module docstring

    bind_simulation<transmission_time_gamma>(handle);
    bind_simulation<transmission_time_weibull>(handle);
    bind_simulation<transmission_time_lognormal>(handle);
    bind_simulation<transmission_time_exponential>(handle);

    bind_simulation_lognormal<transmission_time_gamma>(handle);
    bind_simulation_lognormal<transmission_time_weibull>(handle);
    bind_simulation_lognormal<transmission_time_lognormal>(handle);
    bind_simulation_lognormal<transmission_time_exponential>(handle);

    bind_simulation_powerlaw<transmission_time_gamma>(handle);
    bind_simulation_powerlaw<transmission_time_weibull>(handle);
    bind_simulation_powerlaw<transmission_time_lognormal>(handle);
    bind_simulation_powerlaw<transmission_time_exponential>(handle);

    handle.def("depleted_distribution",&depleted_distribution,
        py::arg("SIZE"),
        py::arg("SIM"),
        py::arg("seed"),
        ""
    );

    handle.def("simulation_discrete_zk",&simulation_discrete_zk,
        py::arg("graph"),
        py::arg("SIZE"),
        py::arg("SIM"),
        py::arg("seed"),
        "returns a vector of vectors (N x kmax) of the number of newly infected with degree k at each discrete step on a BA network"
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

    handle.def("simulate_average", &run_simulation_average, 
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

        "Compute the AVERAGE trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
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




    // py::class_<transmission_time>(handle, "TransmissionTime")
    //     .def(py::init<double>(), py::arg("pinf")
    //     );

    // // Define custom exception classes
    // py::register_exception<std::range_error>(handle, "RangeError");
    // py::register_exception<std::runtime_error>(handle, "RuntimeError");
    // py::register_exception<std::logic_error>(handle, "LogicError");

    py::class_<transmission_time_gamma>(handle, "time_distribution_gamma")
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

    py::class_<transmission_time_deterministic>(handle, "time_distribution_deterministic")
        .def(py::init<double>(), py::arg("tau"),
        "\n"
        "Args:\n"
        "   deterministic time tau.\n"
        "\n"
        )
        .def_readonly("tau", &transmission_time_deterministic::value);


    py::class_<transmission_time_weibull>(handle, "time_distribution_weibull")
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

    py::class_<transmission_time_lognormal>(handle, "time_distribution_lognormal")
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
