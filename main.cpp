#include <iostream>
#include <pybind11/pybind11.h>
#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"


namespace py = pybind11;
rng_t engine;


double f(int r0,int size){
    
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool EDGES_CONCURRENT = true;

    std::uniform_int_distribution<> dis(0, size-1);
    std::uniform_real_distribution<> ran(0.0,1.0);

    erdos_reyni network(size,r0,engine);

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,false);

    for (node_t i = 0; i < 1; i++)
    {
        node_t rand_node = dis(engine);
        // const double thermal =-log(1-ran(engine))/Lambda;
        // simulation.add_thermal_infections({ std::make_pair(rand_node, 0)},Lambda_r2,engine);
        simulation.add_infections({ std::make_pair(rand_node, 0)});
    }
                        
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;
        all_times.push_back(point->time);
        
    }
    
}




PYBIND11_MODULE(episimpy, handle) {
    handle.doc() = "pybind11 example plugin"; // optional module docstring

    handle.def("f", &f, "A function that adds two numbers");
}