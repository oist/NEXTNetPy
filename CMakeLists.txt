cmake_minimum_required(VERSION 3.16)
project(episimpy VERSION 0.1.0)

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
    extern/epidemics/epidemics/utility.cpp)

# include_directories(extern/epidemics/epidemics)
include_directories(AFTER SYSTEM extern/epidemics/epidemics extern/epidemics/ext/dyndist)

if(UNIX AND NOT APPLE)
    # Code to execute on Linux
    message("Including pybind11 as a submodule")
    add_subdirectory(extern/pybind11)
    pybind11_add_module(episimpy ${LIB_SOURCES} main.cpp)
    find_package(Boost REQUIRED)
    link_libraries(Boost::boost)
elseif(APPLE)
    find_package(Python COMPONENTS Interpreter Development)
    find_package(pybind11 CONFIG)
    # pybind11 method:
    pybind11_add_module(episimpy
    main.cpp
    ${SIMULATOR_SOURCES}
    ${LIB_SOURCES}
    ${SIMULATOR_HEADERS}
    ${LIB_HEADERS}
    ${LIB_PCH}
    )
endif()