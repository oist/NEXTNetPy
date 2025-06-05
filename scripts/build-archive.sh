#!/bin/bash
set -e

archive=$1
if [ "$archive" == "" ]; then
	echo "Usage: $0 archive" >&2
	exit 1
fi

tar -czf "$archive" \
	README.md \
	COPYING \
	pyproject.toml \
	CMakeLists.txt \
	src/* \
	extern/pybind11/CMakeLists.txt \
	extern/pybind11/tools/*.cmake \
	extern/pybind11/include/* \
	extern/NEXTNet/nextnet/* \
	extern/NEXTNet/nextnet/pstream/* \
	extern/NEXTNet/ext/dyndist/dyndist/* \
	extern/NEXTNet/ext/prio_queue/prio_queue.hpp \
	extern/boost-*/include/*
