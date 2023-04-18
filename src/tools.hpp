#pragma once 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"

namespace py=pybind11;

double euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi);