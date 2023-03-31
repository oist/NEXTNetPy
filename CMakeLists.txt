cmake_minimum_required(VERSION 3.0.0)
project(episimpy VERSION 0.1.0)

add_subdirectory(extern/pybind11)
add_subdirectory(extern/epidemics)
pybind11_add_module(episimpy main.cpp)

# # Boost headers are mandatory, iostreamd and filesystem are optional
# # and allow plots to be generated during test runs
# find_package(Boost REQUIRED OPTIONAL_COMPONENTS iostreams filesystem)
# link_libraries(Boost::boost)
# if (Boost_FILESYSTEM_FOUND)
#     message("Boost filesystem: Found")
# else()
#     message("Boost filesystem: Not Found")
# endif()
# if (Boost_IOSTREAMS_FOUND)
#     message("Boost iostreams: Found")
# else()
#     message("Boost iostreams: Not Found")
# endif()