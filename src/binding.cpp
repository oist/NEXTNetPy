#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/NextReaction.h"
#include "nextnet/temporal_network.h"
#include "nextnet/weighted_network.h"

#include "networkx.hpp"

namespace py = pybind11;

PYBIND11_MODULE(nextnet, handle) {

    handle.doc() = "nextnet module to efficiently simulate an epidemic on any networkx network, custom network, or temporal network.";

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

    //---------------------------------
    //--------SIMULATION OBJECT--------
    //-------FOR STATIC NETWORKS-------
    //---------------------------------

    py::class_<simulate_next_reaction>(handle,"simulation")
        .def(py::init([](network& nw,
                        transmission_time& psi,
                        transmission_time* rho,
                        bool SIR) {

                simulate_next_reaction::params p;
                p.shuffle_neighbours = false;
                p.edges_concurrent = true;
                p.SIR = SIR;
                simulate_next_reaction sim(nw, psi, rho, p);
                return sim;
            }),
            py::arg("network"),
            py::arg("infection_time"),
            py::arg("recovery_time")=nullptr,
            py::arg("SIR")=true,
            py::keep_alive<1,2>(),
            py::keep_alive<1,3>(),
            py::keep_alive<1,4>(),
            ""
        )

        .def("add_infections",
            [](simulate_next_reaction &self,
                const std::vector<std::pair<int,double>>& infections)
             {
                 // you could validate or transform here
                 self.add_infections(infections);
             },
            py::arg("infections_list"),
            "Add a batch of initial or outside infections; pass a list of (node, time) pairs.")

        .def("run",
                [](simulate_next_reaction &self, py::dict opts, rng_t &engine) {
                    // pull out options or use defaults
                    double max_time = opts.contains("time")
                        ? opts["time"].cast<double>()
                        : std::numeric_limits<double>::infinity();

                    int max_steps = opts.contains("total_infected")
                        ? opts["total_infected"].cast<int>()
                        : std::numeric_limits<int>::max();

                    int threshold = opts.contains("max_infected")
                        ? opts["max_infected"].cast<int>()
                        : static_cast<int>(1e10);

                    // containers
                    std::vector<std::tuple<double,int,int,int>> trajectory;
                    std::vector<double> times;
                    std::vector<int> infected_traj;
                    std::vector<int> recovered_traj;

                    // state counters
                    int current_nb_infected = 0;
                    int cumul_nb_infected  = 0;
                    int current_nb_recovered = 0;

                    // main loop
                    while (true) {
                        auto point = self.step(engine);
                        if (!point)
                            break;

                        int event_type = 0;
                        switch (point->kind) {
                            case epidemic_event_kind::outside_infection:
                            case epidemic_event_kind::infection:
                                current_nb_infected++;
                                cumul_nb_infected++;
                                break;
                            case epidemic_event_kind::reset:
                                // a “recovery” event
                                current_nb_infected--;
                                current_nb_recovered++;
                                event_type = 1;
                                break;
                            default:
                                break;
                        }

                        // stop criteria
                        if (point->time   > max_time ||
                            cumul_nb_infected > max_steps ||
                            current_nb_infected > threshold)
                        {
                            break;
                        }

                        // record
                        times.push_back(point->time);
                        infected_traj.push_back(current_nb_infected);
                        recovered_traj.push_back(current_nb_recovered);
                        trajectory.emplace_back(
                            point->time,
                            point->node,
                            point->source_node,
                            event_type
                        );
                    }

                    // pack into a Python dict
                    py::dict result;
                    result["time"]      = times;
                    result["infected"]  = infected_traj;
                    result["recovered"] = recovered_traj;
                    result["data"]      = trajectory;
                    return result;
                },
            py::arg("options"),
            py::arg("engine"),
            R"doc(
                Simulate the epidemic process until a stopping condition is reached.

                Parameters:
                options (dict): Optional control parameters:
                    - "time" (float): maximum simulation time (default: ∞)
                    - "total_infected" (int): stop after this many cumulative infections (default: ∞)
                    - "max_infected" (int): stop if the number of current infections exceeds this (default: ∞)
                engine (rng_t): Random number generator

                Returns:
                list of tuples (time, node, source_node, event_type)
            )doc"
    );
            

    //---------------------------------
    //------------NETWORKS-------------
    //---------------------------------


    py::class_<network>(handle, "network", py::multiple_inheritance())
        .def("nodes", &network::nodes)
        .def("neighbour", &network::neighbour)
        .def("outdegree", &network::outdegree)
        .def("is_undirected", &network::is_undirected)
        .def("is_simple", &network::is_simple)
        .def("adjacencylist",[](network &self){
            const int n = self.nodes();
            std::vector<std::vector<int>> adj;
            adj.reserve(n);
            for (int node = 0; node<n ; node++){
                std::vector<int> neigh;
                for (int idx = 0;; idx++){
                    const int v =self.neighbour(node,idx);
                    if (v == -1)
                        break;
                    neigh.push_back(v);
                }
                adj.push_back(neigh);
            }
            return adj;
        });

    py::class_<weighted_network, network>(handle, "weighted_network", py::multiple_inheritance())
        .def("weighted_adjacencylist",
            [](weighted_network &self) {
                const int n = self.nodes();
                std::vector<std::vector<std::pair<int,double>>> adj;
                adj.reserve(n);
                for (int node = 0; node < n; ++node) {
                    std::vector<std::pair<int,double>> neigh;
                    for (int idx = 0; ; ++idx) {
                        double weight = 0.0;
                        int v = self.neighbour(node, idx, &weight);
                        if (v == -1)
                            break;
                        neigh.emplace_back(v, weight);
                    }
                    adj.emplace_back(std::move(neigh));
                }
                return adj;
            },
            R"pbdoc(
                weighted_adjacencylist() -> List[List[Tuple[int, float]]]

                Compute and return the weighted adjacency list of this network.

                Returns:
                    A list `adj` of length `n` (number of nodes), where each
                    `adj[i]` is itself a list of `(neighbor, weight)` pairs:
                    
                    - `neighbour` (int): index of a node connected to `i`
                    - `weight`   (float): weight of the edge (i → neighbour)

                Example:
                    >>> wnet = weighted_network(...)
                    >>> adj = wnet.weighted_adjacencylist()
                    >>> for nbr, w in adj[0]:
                    ...     print(f"0 → {nbr} with weight {w}")
            )pbdoc"
        );

    py::class_<networkx,network>(handle,"networkx",py::multiple_inheritance())
        .def(py::init<py::list>())
        .def(py::init<py::object>(),py::arg("networkx_graph"),
        ""
    );
    //---------------------------------
    //---------NETWORK MODELS----------
    //---------------------------------
    py::class_<watts_strogatz,network>(handle, "watts_strogatz", py::multiple_inheritance())
        .def(py::init<node_t, int, double, rng_t&>(), py::arg("size"), py::arg("k"), py::arg("p"), py::arg("rng"));

    py::class_<erdos_renyi,network>(handle, "erdos_renyi", py::multiple_inheritance())
        .def(py::init<int, double, rng_t&>(), py::arg("size"), py::arg("average_degree"), py::arg("rng"));

    py::class_<barabasi_albert,network>(handle, "barabasi_albert", py::multiple_inheritance())
        .def(py::init<int, rng_t&, int>(), py::arg("size"), py::arg("rng"), py::arg("m")=1);

    //---------------------------------
    //--------SIMULATION OBJECT--------
    //------FOR TEMPORAL NETWORKS------
    //---------------------------------

    py::class_<simulate_on_temporal_network>(handle, "simulation_temporal")
        .def(py::init([](temporal_network& nw,
                        transmission_time& psi,
                        transmission_time* rho,
                        bool SIR) {
                simulate_next_reaction::params p;
                p.shuffle_neighbours = false;
                p.edges_concurrent  = true;
                p.SIR = SIR;
                auto nr = std::make_unique<simulate_next_reaction>(nw, psi, rho, p);
                // note: we heap‐allocate nr so its address stays valid
                return simulate_on_temporal_network(*nr.release());
            }),
            py::arg("temporal_network"),
            py::arg("infection_time"),
            py::arg("recovery_time") = nullptr,
            py::arg("SIR")= true,
            // keep python engine alive as long as the wrapper exists:
            py::keep_alive<1, 2>(),
            py::keep_alive<1, 3>(),
            py::keep_alive<1, 4>()
        )
        
        .def("add_infections",
            [](simulate_on_temporal_network &self,
                const std::vector<std::pair<int,double>>& infections)
             {
                 // you could validate or transform here
                 self.simulation.add_infections(infections);
             },
            py::arg("infections_list"),
            "Add a batch of initial or outside infections; pass a list of (node, time) pairs.")

        .def("run",
                [](simulate_on_temporal_network &self, py::dict opts, rng_t &engine) {
                    // pull out options or use defaults
                    double max_time = opts.contains("time")
                        ? opts["time"].cast<double>()
                        : std::numeric_limits<double>::infinity();

                    int max_steps = opts.contains("total_infected")
                        ? opts["total_infected"].cast<int>()
                        : std::numeric_limits<int>::max();

                    int threshold = opts.contains("max_infected")
                        ? opts["max_infected"].cast<int>()
                        : static_cast<int>(1e10);
                    
                    bool network_events  = opts.contains("network_events")
                        ? opts["network_events"].cast<bool>()
                        : false;

                    bool epidemic_events = opts.contains("epidemic_events")
                        ? opts["epidemic_events"].cast<bool>()
                        : true;

                    // containers
                    std::vector<std::tuple<double,int,int,int>> trajectory;
                    std::vector<double> times;
                    std::vector<int> infected_traj;
                    std::vector<int> recovered_traj;

                    // state counters
                    int current_nb_infected = 0;
                    int cumul_nb_infected  = 0;
                    int current_nb_recovered = 0;

                    // main loop
                    while (true) {

                        // std::optional<network_or_epidemic_event_t> any_ev = env.simulator -> step(engine,TMAX);
                        std::optional<network_or_epidemic_event_t> any_ev = self.step(engine, max_time);

                        if (any_ev.has_value()) {
                            if (std::holds_alternative<epidemic_event_t>(*any_ev)) {
                                /* Epidemic event */
                                const auto& ev = std::get<epidemic_event_t>(*any_ev);
                                int event_type = 0;
                                switch (ev.kind) {
                                    case epidemic_event_kind::outside_infection:
                                    case epidemic_event_kind::infection:
                                        current_nb_infected++;
                                        cumul_nb_infected++;
                                        break;
                                    case epidemic_event_kind::reset:
                                        // a “recovery” event
                                        current_nb_infected--;
                                        current_nb_recovered++;
                                        event_type = 1;
                                        break;
                                    default:
                                        break;
                                }

                                // stop criteria
                                if (ev.time   > max_time || cumul_nb_infected > max_steps || current_nb_infected > threshold)
                                    break;

                                if (epidemic_events){
                                    // record
                                    times.push_back(ev.time);
                                    infected_traj.push_back(current_nb_infected);
                                    recovered_traj.push_back(current_nb_recovered);
                                    trajectory.emplace_back(
                                        ev.time,
                                        ev.node,
                                        ev.source_node,
                                        event_type
                                    );
                                }                         

                            } else if (std::holds_alternative<network_event_t>(*any_ev)) {
                                /* Network event */
                                const auto& ev = std::get<network_event_t>(*any_ev);

                                // stop criteria
                                if (ev.time   > max_time)
                                    break;

                                if (network_events){
                                    // record
                                    int event_type = static_cast<int>(ev.kind) + 1;
                                    
                                    trajectory.emplace_back(
                                        ev.time,
                                        ev.target_node,
                                        ev.source_node,
                                        event_type
                                    );    
                                }

                            } else {
                                throw std::logic_error("unknown event type");
                            }

                        } else { // no events
                            break;
                        }

                    }                    

                    // pack into a Python dict
                    py::dict result;
                    result["time"]      = times;
                    result["infected"]  = infected_traj;
                    result["recovered"] = recovered_traj;
                    result["data"]      = trajectory;
                    return result;
                },
            py::arg("options"),
            py::arg("engine"),
            R"doc(
            Simulate an epidemic spreading over a temporal network using the next-reaction method.

            The simulation runs until one of the following stopping conditions is met:
            - the maximum simulation time is reached,
            - the total number of cumulative infections exceeds a limit,
            - the current number of infected individuals exceeds a threshold.

            Parameters:
            -----------
            options : dict (optional)
                Dictionary of control parameters:
                - "time" (float): Maximum simulation time (default: infinity).
                - "total_infected" (int): Stop the simulation after this many cumulative infections (default: unlimited).
                - "max_infected" (int): Stop if the current number of infected nodes exceeds this threshold (default: unlimited).
                - "network_events" (bool): Whether to record network dynamics (default: False).
                - "epidemic_events" (bool): Whether to record epidemic events (default: True).

            engine : rng_t
                A random number generator instance for stochastic simulation.

            Returns:
            --------
            dict
                A dictionary with the following keys:
                - "time": List of times at which epidemic events occurred.
                - "infected": List tracking the number of currently infected nodes over time.
                - "recovered": List tracking the number of recovered nodes over time (if SIR model).
                - "data": List of tuples (time, node, source_node, event_type), where event_type is:
                    - 0 for infection,
                    - 1 for recovery,
                    - 2 for link addition (network event),
                    - 3 for link removal (network event)
                    - 4 for instantaneous link (network event).
            )doc"

    );


    //---------------------------------
    //-------TEMPORAL NETWORKS---------
    //---------------------------------
    py::class_<temporal_network,network>(handle, "temporal_network", py::multiple_inheritance())
        .def("next", &temporal_network::next)
        .def("step", &temporal_network::step)
        .def("notify_epidemic_event", &temporal_network::notify_epidemic_event);


    py::class_<activity_driven_network,temporal_network>(handle, "activity_driven_network", py::multiple_inheritance())
        .def(py::init<std::vector<double>, double, double, double, rng_t&>(), py::arg("activity_rates"), py::arg("eta"), py::arg("m"), py::arg("recovery_rate"), py::arg("rng"));

    // py::class_<temporal_sirx_network,temporal_network>(handle, "temporal_sirx_network", py::multiple_inheritance())
    //     .def(py::init<network&, double, double, rng_t&>(), py::arg("network"), py::arg("kappa0"), py::arg("kappa"), py::arg("rng"));

    // py::class_<empirical_contact_network, temporal_network>(handle, "empirical_contact_network", py::multiple_inheritance())
    //     .def(py::init([](std::string path_to_file, bool is_finite_duration, double dt) {
    //         empirical_contact_network::edge_duration_kind contact_type = 
    //             is_finite_duration ? empirical_contact_network::finite_duration : empirical_contact_network::infitesimal_duration;
    //         std::fstream file(path_to_file);
    //         if (!file)
    //             throw std::runtime_error("failed to open file" + path_to_file);
    //         return new empirical_contact_network(file, contact_type, dt);
    //     }));


    //---------------------------------
    //---------DISTRIBUTIONS-----------
    //---------------------------------
    
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
    #if 0
    //---------------------------------
    //-----------TOOLS ----------------
    //---------------------------------

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
#endif

}
