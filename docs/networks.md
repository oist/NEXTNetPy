# Networks
This document describes the available network classes for building graphs in NextNet.

## Static networks

All static networks inherit the follow common methods:

* **`nodes() -> int`**
  Number of nodes in the network.

* **`neighbour(node: int, idx: int) -> int`**
  Get the `idx`th neighbor of `node`, or `-1` if none.

* **`outdegree(node: int) -> int`**
  Number of outgoing edges from `node`.

* **`is_undirected() -> bool`**
  True if every edge is bidirectional.

* **`is_simple() -> bool`**
  True if the network has no self‑loops or parallel edges.

* **`adjacencylist() -> List[List[int]]`**
  Returns the full adjacency list: for each node, list of neighbor indices.

* (only if weighted) **`weighted_adjacencylist() -> List[List[Tuple[int, float]]]`**
  Compute the weighted adjacency list. Returns, for each node, a list of `(neighbor_index, weight)` pairs.


### Using networks from the NetworkX library
the NextNet python library  is made so that it operates with the [networkx library](https://networkx.org/) which offers a great way to generate and analyse a large collection of random networks. 

Constructors:

* **`networkx(adjacency_list: list)`**
* **`networkx(nx_graph: networkx.Graph)`**

Wraps a Python NetworkX graph for use in simulations. It works with any networkx graph object, as long as the nodes are labeled as integers.

**Example:**
```python
import networkx as nx
n=10**5
p = 3/n
gx = nx.fast_gnp_random_graph(n,p)
gx=nx.convert_node_labels_to_integers(gx)
g = nn.networkx(gx)
```

---
### Networks from a file
#### `empirical_network`

```python
empirical_network(path: str,
                  undirected: bool = True,
                  simplify: bool = False,
                  idxbase: int = 1,
                  sep: char = ' ')
```

Construct from an edge‐list file.

* **`path`**: path to whitespace‐separated edge list
* **`undirected`**: add edges in both directions (default `True`)
* **`simplify`**: remove self‑loops and parallel edges (default `False`)
* **`idxbase`**: node index base in file (e.g. `0` or `1`, default `1`)
* **`sep`**: field separator (default `' '`)

---

### Random Graph Models

#### `watts_strogatz`

```python
watts_strogatz(size: int,
               k: int,
               p: float,
               rng: rng)
```

Generate a Watts–Strogatz small‑world network.

* **`size`**: number of nodes
* **`k`**: each node connects to `k` nearest neighbors
* **`p`**: rewiring probability
* **`rng`**: random number generator

---

#### `erdos_renyi`

```python
erdos_renyi(size: int,
            average_degree: float,
            rng: rng)
```

Generate an Erdős–Rényi random graph G(n, p) where `p = average_degree/(n-1)`.

* **`size`**: number of nodes `n`
* **`average_degree`**: desired mean degree
* **`rng`**: random number generator

---

#### `barabasi_albert`


```python
barabasi_albert(size: int,
                 rng: rng,
                 m: int = 1)
```

Generate a Barabási–Albert scale‑free network.

* **`size`**: initial number of nodes
* **`rng`**: random number generator
* **`m`**: new edges added per added node (default `1`)

---

#### `configuration_model`

```cpp
py::class_<configuration_model, network>(handle, "configuration_model", py::multiple_inheritance())
```

```python
configuration_model(degreelist: List[int],
                    rng: rng)
```

Generate a random network with a prescribed degree sequence.

* **`degreelist`**: list of desired node degrees
* **`rng`**: random number generator

**Remark:** The sum of the degree list needs to be even.

---

#### `configuration_model_clustered`

This allows the generated networks from a prescribed degree sequence and tune the clustering. This is a personal implementation of the algorithm presented in the article by [Serrano and Boguñá: Tuning clustering in random networks with arbitrary degree distributions](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.036133). If you use this function in your research, please cite the original paper.

Three constructors:

1. **Explicit triangles per class**:

   ```python
   configuration_model_clustered(
       degrees: List[int],
       triangles: List[int],
       beta: float,
       engine: rng)
   ```

2. **Clustering function c(k)**:

   ```python
   configuration_model_clustered(
       degrees: List[int],
       ck: Callable[[int], float],
       beta: float,
       engine: rng)
   ```

3. **Alpha exponent form**:

   ```python
   configuration_model_clustered(
       degrees: List[int],
       alpha: float,
       beta: float,
       engine: rng)
   ```

* **`triangles_unsatisfied`** (`bool`): indicates if requested triangles count could not be met.

* See details in original paper.

## Temporal networks

All temporal networks inherit from the base `temporal_network` interface and extend the static `network` functionality.

---

#### `activity_driven_network`

A mechanistic temporal network model where nodes activate stochastically and form temporary links.

```python
activity_driven_network(
    activity_rates: List[float],
    eta: float,
    m: float,
    recovery_rate: float,
    rng: rng
)
```

* **`activity_rates`** (`List[float]`): Activation rate for each node.
* **`eta`** (`float`): Scaling factor for link‐formation probability.
* **`m`** (`float`): Number of links created per activation event.
* **`recovery_rate`** (`float`): Recovery rate for SIR/SIS dynamics when used as a coupled network.
* **`rng`** (`rng`): Random‐number generator for stochastic activations.

---

#### `empirical_temporal_network`

Constructs a temporal network from empirical contact‐sequence data (edge list with timestamps).

```python
empirical_temporal_network(
    file: str,
    finite_duration: bool,
    dt: float
)
```

* **`file`** (`str`): Path to the contact‐list file. Each line should specify a contact event (e.g., `time, node_i, node_j`).
* **`finite_duration`** (`bool`): If `True`, treat each listed contact as lasting until its next appearance; if `False`, treat contacts as instantaneous at resolution `dt`.
* **`dt`** (`float`): Time resolution for discretizing instantaneous contacts.

---

#### `brownian_proximity_network`

*Not available on the Python library yet (v.1.0). It is accessible on the R library.*

#### `temporal_sirx_network`

*Not available on the Python library yet (v.1.0). It is accessible on the R library.*