#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"

namespace py = pybind11;
rng_t engine;

// Wrapper to convert a python networkx graph to our c++ graph object.
class networkx : public virtual graph_adjacencylist {
public:
    //Initialise graph with an adjacency list
    networkx(py::list py_list){

    for (auto row : py_list) {
        std::vector<int> cpp_row;
        for (auto element : row) {
            cpp_row.push_back(py::cast<node_t>(element));
        }
        adjacencylist.push_back(cpp_row);
    }

    }

    // OR Initialise graph with networkx graph
    networkx(py::object graph){

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
};


double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi){
// double euler_lotka_growth_rate(py::object graph,double mean, double variance){
    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    double r = assortativity(network);

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;

    for (node_t node = 0; node < SIZE; node++ ){
        double k = (double) network.outdegree(node);
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3) ;
    }
    const double MU = ( (1-r) * (k2/k1 - 1 ) + r *(k3/k2 - 1) );
    const double SHAPE = psi.mean * psi.mean / psi.variance;
    const double SCALE = psi.variance / psi.mean;
    // const double SHAPE = mean * mean / variance;
    // const double SCALE = variance / variance;

    //Only valid for Gamma distribution 
    const double GROWTH_RATE = 1/SCALE * ( pow(MU,1/SHAPE) - 1 );

    return GROWTH_RATE;
}

// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<int>> run_simulation(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= true,int INITIAL_INFECTED=1, int seed = 1){
    
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);

    std::vector<double> time_trajectory({});
    std::vector<int> infected_trajectory({});

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);


    for (node_t i = 0; i < INITIAL_INFECTED; i++)
    {
        node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
    }
    
    int current_number_infected = INITIAL_INFECTED;
                        
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;

        switch (point->kind) {
            case event_kind::outside_infection:
            case event_kind::infection:
                current_number_infected++;
                break;
            case event_kind::reset:
                current_number_infected--;
                break;
            default:
                break;
        }
        if (current_number_infected<=0 || point->time > TMAX)
            break;
        time_trajectory.push_back(point->time);
        infected_trajectory.push_back(current_number_infected);
    }
    return std::make_tuple(time_trajectory, infected_trajectory);
}



// Wrapper to run a simulation given a network and transmission distribution
std::tuple<std::vector<double>, std::vector<double>> run_simulation_average(py::object graph,transmission_time_gamma psi, transmission_time_gamma* rho= nullptr,bool SIR=false,double TMAX = 1000, bool EDGES_CONCURRENT= false,int INITIAL_INFECTED=1, int seed = 1, int NB_SIMULATIONS = 1, bool TRIM = false){
    
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    bool SHUFFLE_NEIGHBOURS;
    if (EDGES_CONCURRENT || SIR)
        SHUFFLE_NEIGHBOURS=false;
    else 
        SHUFFLE_NEIGHBOURS = true;

    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);

    // Make a vector of pair, which is more convienient once we need to sort the times.
    std::vector<std::pair<double,double>> time_trajectory;
    std::vector<double> infected_trajectory;

    for (int sim = 0; sim < NB_SIMULATIONS; sim++)
        {
            // std::cout << sim << "/" << NB_SIMULATIONS << "\r";
            py::print(sim,"/",NB_SIMULATIONS,"\r",py::arg("end") = "");
            // Initialise/reset simulation object
            simulate_next_reaction simulation(network, psi,rho,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,SIR);
            
            for (node_t i = 0; i < INITIAL_INFECTED; i++)
            {
                const node_t random_node = uniform_node_distribution(engine);
                // sample without replacement:
                if (simulation.is_infected(random_node)){
                    i--;
                    continue;
                }
                simulation.add_infections({ std::make_pair(random_node, 0)});
            }
            
            double infection_counter = INITIAL_INFECTED / NB_SIMULATIONS;
                                
            for (int step = 0; step < SIZE*NB_SIMULATIONS;step++) {
                auto point = simulation.step(engine);
                
                if (!point )
                    break;

                switch (point->kind) {
                    case event_kind::outside_infection:
                    case event_kind::infection:
                        infection_counter =  1.0 ;
                        break;
                    case event_kind::reset:
                        infection_counter =  - 1.0 ;
                        break;
                    default:
                        break;
                }
                if (point->time > TMAX)
                    break;
                time_trajectory.push_back(std::make_pair(point->time , infection_counter ) );
                // infected_trajectory.push_back(current_number_infected/NB_SIMULATIONS);
            }
        }


    std::vector<double> time_trajectory_trimmed;
    std::vector<double> infected_trajectory_trimmed;

    const int COUNTER = (NB_SIMULATIONS > 1 && TRIM) ? NB_SIMULATIONS : 1;
    const int NORM = (NB_SIMULATIONS > 1 && TRIM) ? 1 : NB_SIMULATIONS;

    sort(time_trajectory.begin(), time_trajectory.end(),[](const auto& t1, const auto& t2) {
        return t1.first < t2.first;
    });

    for (int index = 0; index < time_trajectory.size(); index += COUNTER){
            time_trajectory_trimmed.push_back(time_trajectory[index].first);
            infected_trajectory_trimmed.push_back(time_trajectory[index].second / NORM);
    }

    // compute the cumulative sum
    std::partial_sum(infected_trajectory_trimmed.begin(), infected_trajectory_trimmed.end(), infected_trajectory_trimmed.begin());

    return std::make_tuple(time_trajectory_trimmed, infected_trajectory_trimmed);
}






PYBIND11_MODULE(episimpy, handle) {
    handle.doc() = "episimpy module to efficiently simulate an epidemic on any networkx graph."; // optional module docstring

    handle.def("growth_rate",&euler_lotka_growth_rate,
        py::arg("graph"),
        py::arg("infection_time"),
        "Computes the growth rate that solves the generalised Euler-Lotka equation for a SI epidemic with a gamma distribution.\n"
        "\n"
        "Args:\n"
        "    graph (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (episimpy.time_distribution): a gamma distribution object that describes the infection times.\n"
        "\n"
        "Returns:\n"
        "    growth_rate (float): the predicted growth rate of the epidemic.");


    handle.def("simulate", &run_simulation, 
        py::arg("graph"),
        py::arg("infection_time"),
        py::arg("recovery_time")=nullptr,
        py::arg("SIR")=false,
        py::arg("TMAX")=1000,
        py::arg("concurrent_edges")=true,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        "Compute the AVERAGE trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
        "\n"
        "Args:\n"
        "    g (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (episimpy.time_distribution): a distribution object that describes the (attempted) infection times.\n"
        "   recovery_time (episimpy.time_distribution) = None. a dist. obj. that describes the times of recovery. \n"
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
        py::arg("concurrent_edges")=false,
        py::arg("initial_infected")=1,
        py::arg("seed")=0, 
        py::arg("nb_simulations")=1,
        py::arg("trim")=true,

        "Compute the AVERAGE trajectory of a SI/SIR/SIS epidemic on a SINGLE graph G.\n"
        "\n"
        "Args:\n"
        "    g (nx.graph): graph generated by the Networkx module.\n"
        "    infection_time (episimpy.time_distribution): a distribution object that describes the (attempted) infection times.\n"
        "   recovery_time (episimpy.time_distribution) = None. a dist. obj. that describes the times of recovery. \n"
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

    py::class_<transmission_time_gamma>(handle, "time_distribution")
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


}