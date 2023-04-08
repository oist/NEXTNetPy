# episimpy

## Set up

Clone the repository and add the submodules

```
git clone git@github.com:samuelcure/episimpy.git
git submodule update --init --recursive
```

On Linux pybind11 is included as a submodule, however if you are on MacOS you need to install pybind11

```
HOMEBREW_NO_AUTO_UPDATE=1 brew install pybind11
```

Make sure that boost is also installed in your system:

```
HOMEBREW_NO_AUTO_UPDATE=1 brew install boost
```

build the module:

```
mkdir -p build
cd build

cmake ..
make -j4
```

Now if you open a python interpreter inside the build folder you can use the library.
Make sure that networkx is installed before hand:

```
pip install networkx
```

or using anaconda

```
conda install -c anaconda networkx 
```

## How to use

```python
import episimpy
import networkx as nx

n = 10**5 # size of the network
m = 1 # number of edges added per node
G = nx.barabasi_albert(n,m)

# Define the distribution for the infection times
psi = episimpy.time_distribution(10,1)

# Define the distribution for the recovery times
rho = e.time_distribution(15,1)

# To simulate a SI epidemic
times, infected = episimpy.simulate(g,psi)

# To simulate a SIR epidemic
times, infected = episimpy.simulate(g,psi,rho,SIR=True)

# To simulate a SIS epidemic
# WARNING: make sure to define a maximum time since the epidemic
# might never reach an absorbing state.
# by default TMAX=1000, which is VERY long.

times, infected = episimpy.simulate(g,psi,rho,TMAX = 100)

```
