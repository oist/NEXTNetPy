cmake_minimum_required(VERSION 3.0.0)
project(episimpy VERSION 0.1.0)

add_subdirectory(extern/pybind11)
add_subdirectory(extern/epidemics)
pybind11_add_module(episimpy main.cpp)

