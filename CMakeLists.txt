cmake_minimum_required(VERSION 3.5)
project(nmepinet VERSION 0.1.0)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# *** Simulation library sources and headers
set(LIB_SOURCES
    extern/epidemics/epidemics/graph.cpp
    extern/epidemics/epidemics/dynamic_graph.cpp
    extern/epidemics/epidemics/algorithm.cpp
    extern/epidemics/epidemics/random.cpp
    extern/epidemics/epidemics/NextReaction.cpp
    extern/epidemics/epidemics/NextReactionMeanField.cpp
    extern/epidemics/epidemics/nMGA.cpp
    extern/epidemics/epidemics/REGIR.cpp
    extern/epidemics/epidemics/utility.cpp)

set(PY_SOURCES
    src/binding.cpp
    src/networkx.cpp
    src/measure_performance.cpp
    src/simulation_wrapper.cpp
    src/tools.cpp
    src/figures_paper.cpp)


# include_directories(extern/epidemics/epidemics)
include_directories(AFTER SYSTEM extern/epidemics/epidemics extern/epidemics/ext/dyndist)

# Code to execute on Linux
message("Including pybind11 as a submodule")
add_subdirectory(extern/pybind11)
pybind11_add_module(nmepinet ${LIB_SOURCES} ${PY_SOURCES})
find_package(Boost CONFIG REQUIRED)
link_libraries(Boost::boost)
find_package(TBB REQUIRED)
target_link_libraries(nmepinet PRIVATE TBB::tbb)
include_directories(${Boost_INCLUDE_DIRS})
