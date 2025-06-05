# NEXT-Net

*NEXT-Net* is a python module can simulate non-Markovian epidemics on arbitrary networks. *NEXT-Net* is based on a pybind11 wrapper around the C++ simulator [link](https://github.com/oist/NEXTNet)

## Overview
NEXT-Net allows users to efficiently simulate epidemics (SI/SIR/SIS) on networks using various distributions (gamma, Weibull, lognormal, exponential,deterministic).

The main algorithms and their performance are discussed in our [preprint](https://arxiv.org/abs/2412.07095).

## Installation

Download the [latest released](https://github.com/oist/NEXTNetPy/releases) version of *NEXTNetPy-v\<version\>-pkg.tar.gz* and install with

    pip install NEXTNetPy-v<version>-pkg.tar.gz
   
Since *NEXT-Net* is implemented in C++, a C++ compiler is required to install *NEXTNetPy*. Alternatively, if [Git](https://git-scm.com/downloads) is available, the [latest released](https://github.com/oist/NEXTNetPy/releases) version of *NEXTNetPy* can be downloaded, built and installed with

    git clone --recurse-submodules --branch latest-release https://github.com/oist/NEXTNetPy.git
    cd NEXTNetPy
    pip install .    

## Getting Started
To get started with NEXT-Net, check out the [Getting Started](getting_started.md) guide for installation and basic usage instructions.

## Documentation Sections
- [Simulations](simulations.md): Learn how to set up and run simulations.
- [Transmission Distributions](transmission_distributions.md): Overview of available transmission models.
- [Networks](networks.md): Information on network structures supported by NEXT-Net.

## Tutorials & Examples
Explore the example simulations:
- [Basic SI/SIR/SIS Simulations](examples/simulations_examples.md)
- [Temporal Simulations](examples/simulations_temporal.md)

## Contact & Support
For questions, support, or if you think you found a bug, feel free to reach out via GitHub Issues or email samuel.cure@oist.jp


*Note: The documentation is currently being written, typos and small mistakes may be remaining.*

TEST
