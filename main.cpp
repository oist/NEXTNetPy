#include <pybind11/pybind11.h>
// #include "extern/epidemics/epidemics/NextReaction.h"
#include <iostream>


namespace py = pybind11;

double f(double x, double y){
    return x+y;
}




PYBIND11_MODULE(episimpy, handle) {
    handle.doc() = "pybind11 example plugin"; // optional module docstring

    handle.def("f", &f, "A function that adds two numbers");
}