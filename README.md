# episimPY

*episimpy* is a python module that allows non-Markovian epidemics for arbitrary transmission time and recovery (reset) time distribution to be simulated on arbitrary contact network graphs. episimpy is based on a pybind11 wrapper around the C++ simulator https://github.com/samuelcure/Epidemics-On-Networks.

## Set up

Open a new terminal tab, clone the repository and add the submodules

```
git clone git@github.com:samuelcure/episimpy.git
git submodule update --init --recursive
```

On Linux pybind11 is included as a submodule, however if you are on MacOS you may need to install pybind11 globally using homebrew

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

# Define the distribution for the infection times (Gamma distributed)
MEAN_INFECTION = 5
VARIANCE_INFECTION = 1
psi = episimpy.time_distribution(MEAN_INFECTION,VARIANCE_INFECTION)

# Define the distribution for the recovery times
MEAN_RECOVER = 7
VARIANCE_RECOVERY= 1
rho = episimpy.time_distribution(MEAN_RECOVERY,VARIANCE_RECOVERY)

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

## Parameters

the `episimpy.simulate` can take the following parameters:
* `recovery_time` A gamma distribution, by default is `None`
* `SIR` by default is `False`, and has an effect only if the recovery_time is defined
* `concurrent_edges` by default is set to `False`. The current implementation adds edges sequentially in the priority queue, since some edges despite being active, do not contribute to the epidemic (for example in networks with short loops). For SI/SIR model, the sequential mode is believed to always be faster, independently of the graph. SIS not sure since the edges have to be reshuffled everytime.
* `seed` choose the random seed for reproducibility.
* `TMAX` maximum time set by default to 1000. Advised to modify when simulating a SIS epidemic.
* `initial_infected` set by default to 1.
