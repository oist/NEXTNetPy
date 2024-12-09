#!/bin/bash

set -o pipefail
set -e

mkdir /src
cd /src
git clone git@github.com:oist/NEXTNetPy.git NEXTNetPy
cd NEXTNetPy
git submodule update --init --recursive
pip install .

cd /
mkdir -p out/
python ./test.py
ls -l /out
