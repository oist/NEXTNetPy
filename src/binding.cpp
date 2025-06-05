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
        .def(py::init<int>(), py::arg("seed"),
            R"(
            **rng(seed: int = None)**

            A reproducible random‐number generator for stochastic simulations.

            Args:
                seed (int, optional): If provided, initializes the generator to a
                    fixed state for reproducible outcomes. If omitted, uses a
                    stochastic seed.
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
                 self.add_infections(infections);
             },
            py::arg("infections_list"),
            "Add a batch of initial or outside infections; pass a list of (node, time) pairs.")

        .def("run",
                [](simulate_next_reaction &self,rng_t &engine, py::dict opts) {
                    double max_time = opts.contains("time")
                        ? opts["time"].cast<double>()
                        : std::numeric_limits<double>::infinity();

                    int max_steps = opts.contains("total_infected")
                        ? opts["total_infected"].cast<int>()
                        : std::numeric_limits<int>::max();

                    int threshold = opts.contains("max_infected")
                        ? opts["max_infected"].cast<int>()
                        : static_cast<int>(1e10);

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
                                current_nb_infected--;
                                current_nb_recovered++;
                                event_type = 1;
                                break;
                            default:
                                break;
                        }

                        // stop criteria
                        if (point->time >= max_time ||
                            cumul_nb_infected >= max_steps ||
                            current_nb_infected >= threshold)
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
            py::arg("engine"),
            py::arg("options")=py::dict(),
            R"(
            **run(options: dict, engine: rng) -> dict**

            Execute the epidemic simulation until a stopping condition is met.

            **Options keys (all optional):**
            - `"time"` (float): Max simulation time (default ∞)
            - `"total_infected"` (int): Max cumulative infections (default ∞)
            - `"max_infected"` (int): Max concurrent infections (default ∞)

            Args:
                engine (rng): Random‐number generator instance.

            Returns:
                A dict with:
                - `"time"`: List[float] event times  
                - `"infected"`: List[int] current infected counts  
                - `"recovered"`: List[int] recovered counts  
                - `"data"`: List[(time, node, source_node, event_type)]  
            )"
    );
            

    //---------------------------------
    //------------NETWORKS-------------
    //---------------------------------


    py::class_<network>(handle, "network", py::multiple_inheritance())
        .def("nodes",        &network::nodes,
            R"(
            **nodes() -> int**

            Number of nodes in the network.
            )"
        )
        .def("neighbour",    &network::neighbour,
            R"(
            **neighbour(node: int, idx: int) -> int**

            Get the idx-th neighbor of `node`; returns -1 if none.
            )"
        )
        .def("outdegree",    &network::outdegree,
            R"(
            **outdegree(node: int) -> int**

            Number of outgoing edges from `node`.
            )"
        )
        .def("is_undirected",&network::is_undirected,
            R"(
            **is_undirected() -> bool**

            True if every edge is bidirectional.
            )"
        )
        .def("is_simple",    &network::is_simple,
            R"(
            **is_simple() -> bool**

            True if network has no self-loops or parallel edges.
            )"
        )
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
        },
            R"(
            **adjacencylist() -> List[List[int]]**

            Returns the adjacency list: for each node, a list of neighbor indices.
            )"
        );

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
            **weighted_adjacencylist() -> List[List[Tuple[int, float]]**

            Compute the weighted adjacency list.

            Returns:
                A list of length `n` where each entry is a list of
                (neighbor_index, edge_weight) pairs.
            )pbdoc"
        );

    py::class_<networkx,network>(handle,"networkx",py::multiple_inheritance())
        .def(py::init<py::list>())
        .def(py::init<py::object>(),py::arg("networkx_graph"),
        R"(
        **networkx(nx_graph: List or networkx.Graph)**

        Wraps a Python NetworkX graph object for use in simulations.
        )"
    );

    py::class_<empirical_network, network>(handle, "empirical_network", py::multiple_inheritance())
        .def(
            // wrap the std::istream& ctor in a lambda
            py::init([](const std::string& path,
                        bool undirected,
                        bool simplify,
                        int idxbase,
                        char sep) {
                std::ifstream file(path);
                if (!file.is_open()) {
                    throw std::runtime_error("Could not open file: " + path);
                }
                return empirical_network(file, undirected, simplify,
                                         static_cast<node_t>(idxbase),
                                         sep);
            }),
            py::arg("path"),
            py::arg("undirected")  = true,
            py::arg("simplify")    = false,
            py::arg("idxbase")     = 1,
            py::arg("sep")         = ' ',
            R"doc(
                Construct an empirical_network by reading an edge list from a file.

                Parameters
                ----------
                path : str
                    Path to a whitespace‐separated edge list file.
                undirected : bool, optional
                    If true, add edges in both directions.  Default: True.
                simplify : bool, optional
                    If true, remove self‐loops and parallel edges.  Default: False.
                idxbase : int, optional
                    The node index base used in the file (e.g. 1 or 0).  Default: 1.
                sep : char, optional
                    Separator character between fields.  Default: ' '.
            )doc"
        );


    //---------------------------------
    //---------NETWORK MODELS----------
    //---------------------------------

    py::class_<watts_strogatz,network>(handle, "watts_strogatz", py::multiple_inheritance())
        .def(py::init<node_t, int, double, rng_t&>(), py::arg("size"), py::arg("k"), py::arg("p"), py::arg("rng"),
            R"(
            **watts_strogatz(size: int, k: int, p: float, rng: rng)**

            Generate a Watts–Strogatz small‐world network.

            Args:
                size: Number of nodes.
                k: Each node connects to k nearest neighbors.
                p: Rewiring probability.
                rng: RNG generator.
            )"
    );

    py::class_<erdos_renyi,network>(handle, "erdos_renyi", py::multiple_inheritance())
        .def(py::init<int, double, rng_t&>(), py::arg("size"), py::arg("average_degree"), py::arg("rng"),
            R"(
            **erdos_renyi(size: int, average_degree: float, rng: rng)**

            Generate an Erdős–Rényi random graph G(n, p).

            Args:
                size: Number of nodes n.
                average_degree: Desired mean degree (p = average_degree/(n-1)).
                rng: RNG generator.
            )"
    );

    py::class_<barabasi_albert,network>(handle, "barabasi_albert", py::multiple_inheritance())
        .def(py::init<int, rng_t&, int>(), py::arg("size"), py::arg("rng"), py::arg("m")=1,
            R"(
            **barabasi_albert(size: int, rng: rng, m: int = 1)**

            Generate a Barabási–Albert scale‐free network.

            Args:
                size: Initial number of nodes.
                rng: RNG generator.
                m: New edges per added node.
            )"
    );

    py::class_<config_model,network>(handle, "configuration_model", py::multiple_inheritance())
        .def(py::init<std::vector<int>, rng_t&>(), py::arg("degreelist"), py::arg("rng"),
            R"(
            **configuration_model(degreelist: int, rng: rng)**

            Generate a network with a prescribed degree distribution

            Args:
                size: Initial number of nodes.
                rng: RNG generator.
            )"
    );


    py::class_<config_model_clustered_serrano, network>(handle, "configuration_model_clustered",py::multiple_inheritance())
        .def(py::init<
                std::vector<int>,              // degrees
                std::vector<int>,              // triangles
                double,                        // beta
                rng_t&                         // engine
            >(),
            py::arg("degrees"),
            py::arg("triangles"),
            py::arg("beta"),
            py::arg("engine"),
            R"pbdoc(
                Constructor with explicit triangle counts per degree class.

                Args:
                    degrees (List[int]): node degrees.
                    triangles (List[int]): triangles per degree class.
                    beta (float): degree-class probability parameter P(k).
                    engine (rng_t): random number generator.
            )pbdoc"
        )
        .def(py::init<
                std::vector<int>,              // degrees
                std::function<double(int)>,    // ck(k)
                double,                        // beta
                rng_t&                         // engine
            >(),
            py::arg("degrees"),
            py::arg("ck"),
            py::arg("beta"),
            py::arg("engine"),
            R"pbdoc(
                Constructor with clustering function c(k).

                Args:
                    degrees (List[int]): node degrees.
                    ck (Callable[[int], float]): clustering probability function.
                    beta (float): degree-class probability parameter P(k).
                    engine (rng_t): random number generator.
            )pbdoc"
        )
        .def(py::init<
                std::vector<int>,              // degrees
                double,                        // alpha
                double,                        // beta
                rng_t&                         // engine
            >(),
            py::arg("degrees"),
            py::arg("alpha"),
            py::arg("beta"),
            py::arg("engine"),
            R"pbdoc(
                Constructor with alpha exponent for c(k)=0.5*(k-1)^alpha.

                Args:
                    degrees (List[int]): node degrees.
                    alpha (float): exponent parameter for c(k).
                    beta (float): degree-class probability parameter P(k).
                    engine (rng_t): random number generator.
            )pbdoc"
        )

        .def_readwrite(
            "triangles_unsatisfied",
            &config_model_clustered_serrano::triangles_unsatisfied,
            "True if the requested triangles count could not be satisfied."
    );


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
                 self.simulation.add_infections(infections);
             },
            py::arg("infections_list"),
            "Add a batch of initial or outside infections; pass a list of (node, time) pairs.")

        .def("run",
                [](simulate_on_temporal_network &self, rng_t &engine,py::dict opts) {
                    // pull out options or use defaults
                    double max_time = opts.contains("time")
                        ? opts["time"].cast<double>()
                        : std::numeric_limits<double>::infinity();

                    int max_steps = opts.contains("total_infected")
                        ? opts["total_infected"].cast<int>()
                        : static_cast<int>(1e10);

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
                                if ((ev.time >= max_time) ||
                                    (cumul_nb_infected >= max_steps) ||
                                    (current_nb_infected >= threshold))
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
                py::arg("engine"),
                py::arg("options")=py::dict(),
            R"doc(
            **run(options: dict, engine: rng) -> dict**

            Execute the next‐reaction epidemic simulation on a temporal network.

            **Options (all optional):**
            - `"time"` (float): Max simulation time (∞)
            - `"total_infected"` (int): Max cumulative infections (∞)
            - `"max_infected"` (int): Max concurrent infections (∞)
            - `"network_events"` (bool): Record network events? (default False)
            - `"epidemic_events"` (bool): Record epidemic events? (default True)

            Returns a dict:
            - `"time"`: times of epidemic events  
            - `"infected"`: infected count over time  
            - `"recovered"`: recovered count over time  
            - `"data"`: List of (time, node, source_node, event_type) where
            - 0 = infection, 1 = recovery, 2 = link added, 3 = link removed, 4 = instant link.
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
        .def(py::init<std::vector<double>, double, double, double, rng_t&>(), py::arg("activity_rates"), py::arg("eta"), py::arg("m"), py::arg("recovery_rate"), py::arg("rng"),
        R"(
        **activity_driven_network(activity_rates, eta, m, recovery_rate, rng)**

        Temporal network where nodes activate at given rates.

        Args:
            activity_rates (List[float]): Activation rate per node.
            eta (float): Probability scaling for link‐formation.
            m (float): Number of links created per activation.
            recovery_rate (float): Recovery rate (for SIR dynamics).
            rng (rng): RNG for stochastic activations.
        )"    
    );

    // py::class_<temporal_sirx_network,temporal_network>(handle, "temporal_sirx_network", py::multiple_inheritance())
    //     .def(py::init<network&, double, double, rng_t&>(), py::arg("network"), py::arg("kappa0"), py::arg("kappa"), py::arg("rng"));

    py::class_<empirical_contact_network, temporal_network>(handle, "empirical_contact_network", py::multiple_inheritance())
        .def(
            py::init([](const std::string & path_to_file,
                        bool               is_finite_duration,
                        double             dt,
                        double             weight)
            {
                // edge type
                empirical_contact_network::edge_duration_kind contact_type = is_finite_duration
                    ? empirical_contact_network::edge_duration_kind::finite_duration
                    : empirical_contact_network::edge_duration_kind::infitesimal_duration;

                // open file
                auto file = std::make_shared<std::fstream>(path_to_file, std::ios::in);
                if (! file->is_open()) {
                    throw std::runtime_error("failed to open file: " + path_to_file);
                }

                return new empirical_contact_network(*file, contact_type, dt, weight);
            }),
            py::arg("path_to_file"),
            py::arg("finite_duration"),
            py::arg("dt"),
            py::arg("weight") = 1.0,
            R"pbdoc(
                empirical_contact_network(path_to_file: str,
                                         finite_duration: bool,
                                         dt: float,
                                         weight: float = 1.0)

                Build an empirical_contact_network from a contact‐sequence file.

                Parameters
                ----------
                path_to_file : str
                    Path to the file containing contact times (one line per contact).
                finite_duration : bool
                    If True, interpret edges as having a finite duration `dt`; 
                    otherwise, treat them as instantaneous.
                dt : float
                    Time‐step (duration) to use when you specify `finite_duration == True`.
                weight : float, optional (default=1.0)
                    An uniform weight to assign to each edge. Defaults to 1.0.
            )pbdoc"
        )

        .def(
            "compute_number_of_edges",
            &empirical_contact_network::compute_number_of_edges,
            py::arg("engine"),
            R"pbdoc(
                compute_number_of_edges(engine: rng_t) -> List[List[float]]

                Returns the time trajectory of the number of edges present in the network.
            )pbdoc"
        );

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

    py::class_<transmission_time_infectiousness, transmission_time>(handle, "transmission_time_infectiousness")
        .def(py::init<std::vector<double>,std::vector<double>>(), py::arg("tau"),py::arg("lambda"));

    //---------------------------------
    //-----------TOOLS ----------------
    //---------------------------------


    handle.def("assortativity",
        &assortativity,
        py::arg("network"),
        R"pbdoc(
            Compute the assortativity coefficient of a network.

            Args:
                network (network): The network to analyze.

            Returns:
                float: the assortativity coefficient.
        )pbdoc"
    );

    handle.def("reproduction_matrix",
        [](network &nw) {
            double out_r = 0.0, out_c = 0.0;
            double out_k1 = 0.0, out_k2 = 0.0, out_k3 = 0.0;
            double out_m1 = 0.0, out_m2 = 0.0;
            double out_R0 = 0.0, out_R_r = 0.0, R_pert = 0.0;

            auto mat = reproduction_matrix(
                nw,
                &out_r, &out_c,
                &out_k1, &out_k2, &out_k3,
                &out_m1, &out_m2,
                &out_R0, &out_R_r,
                &R_pert
            );

            std::vector<double> params = {
                out_r, out_c,
                out_k1, out_k2, out_k3,
                out_m1, out_m2,
                out_R0, out_R_r,
                R_pert
            };
            return py::make_tuple(mat, params);
        },
        py::arg("network"),
        R"pbdoc(
            Compute the reproduction matrix and return summary parameters.

            Args:
                network (network): The network to analyze.

            Returns:
                tuple:
                - matrix (List[List[float]]): the reproduction matrix.
                - params_list (List[float]): [r, c, k1, k2, k3, m1, m2, R0, R_r, R_pert].
        )pbdoc"
    );


}
