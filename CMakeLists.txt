cmake_minimum_required(VERSION 3.15)
project(nextnet VERSION 0.1.0)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# *** Simulation library sources and headers
set(LIB_SOURCES
    extern/NEXTNet/nextnet/network.cpp
    extern/NEXTNet/nextnet/temporal_network.cpp
	extern/NEXTNet/nextnet/weighted_network.cpp
    extern/NEXTNet/nextnet/algorithm.cpp
    extern/NEXTNet/nextnet/random.cpp
    extern/NEXTNet/nextnet/NextReaction.cpp
    extern/NEXTNet/nextnet/NextReactionMeanField.cpp
    extern/NEXTNet/nextnet/nMGA.cpp
    extern/NEXTNet/nextnet/REGIR.cpp
    extern/NEXTNet/nextnet/utility.cpp)

set(PY_SOURCES
    src/binding.cpp
    src/networkx.cpp)

include_directories(AFTER SYSTEM
	extern/boost-core/include
	extern/boost-config/include
	extern/boost-utility/include
	extern/boost-preprocessor/include
	extern/boost-mpl/include
	extern/boost-static_assert/include
	extern/boost-assert/include
	extern/boost-type_traits/include
	extern/boost-move/include
	extern/boost-throw_exception/include
	extern/boost-tuple/include
	extern/boost-bind/include
	extern/boost-detail/include
	extern/boost-iterator/include
	extern/boost-range/include
	extern/boost-math/include
	extern/boost-container_hash/include
	extern/boost-describe/include
	extern/boost-mp11/include
	extern/boost-lexical_cast/include
	extern/boost-multi_index/include
	extern/NEXTNet
        extern/NEXTNet/ext
	extern/NEXTNet/ext/dyndist)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Found Python version: ${Python3_VERSION}")

# Keep things single-threaded, also avoid problems with std::unary_function
# which was removed with C++ 17
add_compile_definitions(_PSTL_PAR_BACKEND_SERIAL)
add_compile_definitions(NEXTNET_BOOST_NO_CXX17_HDR_EXECUTION)
add_compile_definitions(NEXTNET_BOOST_NO_CXX98_FUNCTION_BASE)

# Disable FORTIFY_SOURCE
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE")

message("Including pybind11 as a submodule")
add_subdirectory(extern/pybind11)
pybind11_add_module(nextnet ${LIB_SOURCES} ${PY_SOURCES})
set_target_properties(nextnet PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/nextnet)
