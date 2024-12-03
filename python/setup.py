from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os

# Define the Boost and other library dependencies
include_dirs = [
    os.path.abspath("../extern/pybind11/include"),  # Pybind11 headers
    os.path.abspath("../extern/epidemics/include"),  # Add custom paths
    os.getenv("BOOST_INCLUDE_DIR", "/usr/include/boost"),  # Boost headers
]

library_dirs = [
    os.getenv("BOOST_LIBRARY_DIR", "/usr/lib"),  # Boost libraries
]

libraries = [
    "boost_system",  # Add necessary Boost libraries
]

ext_modules = [
    Pybind11Extension(
        "nmepinet",  # Replace with your Python module name
        [
            "../src/binding.cpp",  # Entry point for Pybind11 binding
            "../src/networkx.cpp",  # Add other C++ source files
            "../src/tools.cpp",
            "../src/simulation_wrapper.cpp",
        ],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=libraries,
        cxx_std=17,  # Ensure modern C++ standard
    ),
]

setup(
    name="my_project",  # Replace with your desired project name
    version="0.1.0",
    author="Samuel Cure",
    author_email="samuel.cure@oist.jp",
    description="Simulate Epidemics on Complex Networks",
    long_description=open("../README.md").read(),
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
