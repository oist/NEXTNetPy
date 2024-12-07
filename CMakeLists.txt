cmake_minimum_required(VERSION 3.15)
project(nextnet VERSION 0.1.0)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# *** Simulation library sources and headers
set(LIB_SOURCES
    extern/NEXTNet/nextnet/network.cpp
    extern/NEXTNet/nextnet/temporal_network.cpp
    extern/NEXTNet/nextnet/algorithm.cpp
    extern/NEXTNet/nextnet/random.cpp
    extern/NEXTNet/nextnet/NextReaction.cpp
    extern/NEXTNet/nextnet/NextReactionMeanField.cpp
    extern/NEXTNet/nextnet/nMGA.cpp
    extern/NEXTNet/nextnet/REGIR.cpp
    extern/NEXTNet/nextnet/utility.cpp)

set(PY_SOURCES
    src/binding.cpp
    src/networkx.cpp
    src/simulation_wrapper.cpp
    src/tools.cpp)

include_directories(AFTER SYSTEM
	extern/boost-core/include
	extern/boost-math/include
	extern/boost-utility/include
	extern/NEXTNet/nextnet
	extern/NEXTNet/ext/dyndist)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Found Python version: ${Python3_VERSION}")

message("Including pybind11 as a submodule")
add_subdirectory(extern/pybind11)
pybind11_add_module(nextnet ${LIB_SOURCES} ${PY_SOURCES})
set_target_properties(nextnet PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/nextnet)
# TBB needed for me
find_package(TBB REQUIRED)
target_link_libraries(nextnet PRIVATE TBB::tbb)