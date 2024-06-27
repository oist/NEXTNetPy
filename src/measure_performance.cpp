#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "random.h"
#include "analysis.h"
#include "NextReaction.h"
#include "networkx.hpp"
#include "simulation_wrapper.hpp"
#include "tools.hpp"
#include "graph.h"
#include "measure_performance.hpp"


namespace py = pybind11;