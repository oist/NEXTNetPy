# episimpy

## Set up ##
Clone the repository and add the submodules

```
git clone ...
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

## How to use ##

```python
import episimpy
import networkx as nx

n = 10**5 # size of the network
m = 1 # number of edges added per node
G = nx.barabasi_albert(n,m)

episimpy.simulate(n,m)


```
