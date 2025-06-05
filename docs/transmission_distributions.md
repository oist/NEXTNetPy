# Transmission Distributions

List of the current distributions available to generate infection times.

## Available Transmission Distributions

### `Gamma`
A gamma-distributed transmission time model.

**Constructor**
```cpp
transmission_time_gamma(double mean, double variance, double pinf=0.0)
```
- **mean**: The mean transmission time.
- **variance**: The variance of the transmission time.
- **pinf**: (Optional) Probability that the transmission does not occur (default: `0.0`).

**Example:**
```python
psi = nn.transmission_time_gamma(3.0, )
```

### `Weibull `
A Weibull-distributed transmission time model.

**Constructor**
```cpp
transmission_time_weibull(double shape, double scale, double pinf=0.0)
```
- **shape**: The shape parameter of the Weibull distribution.
- **scale**: The scale parameter of the Weibull distribution.
- **pinf**: (Optional) Probability that the transmission does not occur (default: `0.0`).

**Example:**
```python
psi = nn.transmission_time_weibull(2.0, 5.0)
```

### `Lognormal `
A lognormal-distributed transmission time model.

**Constructor**
```cpp
transmission_time_lognormal(double mean, double variance, double pinf=0.0)
```
- **mean**: The mean of the lognormal distribution.
- **variance**: The variance of the lognormal distribution.
- **pinf**: (Optional) Probability that the transmission does not occur (default: `0.0`).

**Example:**
```python
psi = nn.transmission_time_lognormal(3.0, 1.2)
```

### `Exponential`
An exponentially-distributed transmission time model.

**Constructor**
```cpp
transmission_time_exponential(double rate)
```
- **rate**: The rate parameter of the exponential distribution, where mean transmission time is `1/rate`.

**Example:**
```python
psi = nn.transmission_time_exponential(0.5)
```

### `Delta (deterministic)`
A deterministic transmission time model where all transmission times are fixed.

**Constructor**
```cpp
transmission_time_deterministic(double tau)
```
- **tau**: The fixed transmission time.

**Example:**
```python
psi = nn.transmission_time_deterministic(4)
```

### `Infectiousness (custom distribution)`
A transmission time $\psi(\tau)$ that is defined from the infectiousness, or hazard rate $\lambda(\tau)$. The user enters an array `tau` and an array `infectiousness` of same length to represent $\lambda(\tau)$. The resulting distribution is then given by
$$ \psi(\tau)=\lambda(\tau)\exp\left(-\int_0^\tau \lambda(\tau')\mathrm{d}\tau'\right).$$

**Constructor**
```cpp
transmission_time_infectiousness(double tau)
```
- **tau**: The fixed transmission time.

**Example:**
```python
tmax = 10
dt = 0.1
tau = [ dt * j for j in range(int(tmax/dt))]
infectiousness = [3 for t in tau]
psi = nn.transmission_time_infectiousness(tau,infectiousness)
```

**Remark**
- `transmission_time_infectiousness` linearly interpolates the values of the infectiousness.
- Beyond `tmax`, the infectiousness is assumed to be constant and takes the last value of the `infectiousness` array.
- The first value of the `infectiousness` array cannot be zero currently (v.0.4.0).
- This function can be particularly useful for simulations on temporal networks where the user might be more interested in defining the infection times from the infectiousness, or hazard rate, rather than the first infection time attempt.