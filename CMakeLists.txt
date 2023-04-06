cmake_minimum_required(VERSION 3.16)
project(episimpy VERSION 0.1.0)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(extern/pybind11)

# *** Simulation library sources and headers
set(LIB_SOURCES
    extern/epidemics/epidemics/graph.cpp
    extern/epidemics/epidemics/dynamic_graph.cpp
    extern/epidemics/epidemics/algorithm.cpp
    extern/epidemics/epidemics/random.cpp
    extern/epidemics/epidemics/NextReaction.cpp
    extern/epidemics/epidemics/NextReactionMeanField.cpp
    extern/epidemics/epidemics/nMGA.cpp
    extern/epidemics/epidemics/utility.cpp)

# include_directories(extern/epidemics/epidemics)
include_directories(AFTER SYSTEM extern/epidemics/epidemics extern/epidemics/ext/dyndist)
pybind11_add_module(episimpy ${LIB_SOURCES} main.cpp)

find_package(Boost REQUIRED)
link_libraries(Boost::boost)
