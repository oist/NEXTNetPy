#!/bin/bash
tar -czf NEXTNetPy.tar.gz \
	README.md \
	COPYING \
	pyproject.toml \
	CMakeLists.txt \
	src/* \
	extern/pybind11/CMakeLists.txt \
	extern/pybind11/tools/*.cmake \
	extern/pybind11/include/* \
	extern/NEXTNet/nextnet/* \
	extern/NEXTNet/ext/dyndist/dyndist/* \
	extern/NEXTNet/ext/prio_queue/prio_queue.hpp \
	extern/boost-*/include/*
