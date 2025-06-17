# Analysis Tools


---
### `reproduction_matrix`

Compute the reproduction matrix for a given network and return summary parameters. See *Cure, S., Pflug, F. G., & Pigolotti, S. (2025). Exponential rate of epidemic spreading on complex networks. Physical Review E, 111(4), 044311.*

**Constructor**

```python
mat, params_list = nextnet.reproduction_matrix(network)
```

* **network** (`network`): The network to analyze. (needs to be a nextnet object)

**Returns:**

* **mat** (`List[List[float]]`): The reproduction matrix as a nested list. mat[k][k'] is the average number of nodes of degree $k$ a node of degree $k'$ will reach.
* **params\_list** (`List[float]`): A list of 10 floats in order:

  ```
  [r, c, k1, k2, k3, m1, m2, R0, R_r, R_pert]
  ```

  * **r**: *Assortativity: degree correlation of connected nodes.*
  * **c**: *Average number of triangles per node (for nodes of degree $k>1$)*
  * **k1**, **k2**, **k3**: *First three moments of degree distribution*
  * **m1**, **m2**: *$m_1$ is the average number of triangle per edge. $m_2$ is related to the average degree of a node belonging to a triangle, see original paper.*
  * **R0**: *returns $\langle k^2\rangle/\langle k\rangle -1$, even when the network is correlated.*
  * **R\_r**: *Reproduction number including effect of correlation.*
  * **R\_pert**: *Reproduction number including both correlation and clustering effect*

**Example:**

```python
import nextnet as nn
import numpy as np

mat, params = nn.reproduction_matrix(network)
eigenvalues = np.linalg.eigvals(mat)
R = np.max(eigenvalues) # True reproduction number

# Unpack parameters:
r, c, k1, k2, k3, m1, m2, R0, R_r, R_pert = params
print(f"r = {r}, c = {c}, k1 = {k1}, ...")
```

### `generalised_logistic_curve`

A simple implementation of the logistic function using `numpy`. Careful: this function actually returns the log of the generalized logistic function. (currently in v.1.0)


**Function**

```python
def generalised_logistic_curve(t,t0,rate,I0,nu,N):
  """
  The logistic function with three parameters: N, I0, rate.
  
  Parameters:
  - t: Independent variable (array)
  - N: Carrying capacity (size of network)
  - I0: Initial number of infected
  - rate: Exponential growth rate of epidemic
  """

  Q = -1 + np.power(N/I0,nu)
  
  return np.log(N) - 1/nu*np.log(1+Q * np.exp(-rate*nu*(t-t0)))
```



---
### `measure_rate`

Estimate a fitted rate parameter by fitting the log of observed counts against time using a generalized logistic curve. This function is written in python and fits a generalized logistic function to an epidemic trajectory.

**Function**

```python
rate_fit, gen_logi_params = measure_rate(times, infected)
```

* **times** (`array-like`): Sequence of event times (e.g., timestamps). Converted to a NumPy array internally.
* **infected** (`array-like`): Sequence of observed counts corresponding to `times`. Must be positive for log-transform.

**Returns:**

* **rate\_fit** (`float`): Fitted rate computed as `lambda_0 * (1 - eta)`, where `lambda_0` and `eta` derive from the logistic fit.
* **gen\_logi\_params** (`ndarray` of shape `(5,)`): Fitted parameters from `scipy.optimize.curve_fit`, in order:

  ```
  [t0_fit, lambda_0, I0_fit, nu_fit, SIZE_fit]
  ```

  * **t0\_fit**: Time-shift (inflection point) fitted.
  * **lambda\_0**: Fitted exponential growth rate from the logistic model.
  * **I0\_fit**: Fitted initial count.
  * **nu\_fit**: Fitted shape parameter.
  * **SIZE\_fit**: Fitted carrying capacity (initial N guess).

**Example:**

```python
import nextnet as nn
import numpy as np

# simulate epidemic
# times, infected = ....

rate_fit, params = measure_rate(times, infected)
t0_fit, lambda_0, I0_fit, nu_fit, SIZE_fit = params

t = np.array(times)
gen_log_traj = np.exp(nn.generalised_logistic_curve(t,*params))

#plot trajectory
# plt.plot(t,gen_log_traj)

#plot inferred exponential rate
# plt.plot(t,I0_fit * np.exp(t*rate_fit))

```

**Notes & Caveats:**
* This function uses `scipy.optimize.curve_fit` and may issue warnings or errors on unusual or noisy trajectories. For complex cases, consider fitting manually with customized settings.
* Its purpose is to infer the exponential growth rate by fitting the log of the trajectory with equal weighting across time, preventing large counts at later stages from dominating the fit.
* The fitted rate is computed as:

  1. Extract parameters from the logistic fit.
  2. Compute `eta = I0_fit / N_events**nu_fit`.
  3. Return `rate_fit = lambda_0 * (1 - eta)`.
* Because bounds, initial guesses, and weighting may not suit every dataset, adjust or replace this routine when finer control is needed.
