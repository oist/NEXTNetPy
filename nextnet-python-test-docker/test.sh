#!/bin/bash

set -o pipefail
set -e

mkdir /src
cd /src
git clone --branch latest-release https://github.com/oist/NEXTNetPy.git NEXTNetPy
cd NEXTNetPy
git submodule update --init --recursive
pip install .

cd /
mkdir -p out/
rm -f out/test.pdf
python ./test.py
ls -l /out
