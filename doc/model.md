

# Equilibrium Term Structure Models
## Vasicek Model

SDE:

$$
d r_t = \gamma (\bar{r} - r_t) d t + \sigma d X_t
$$
where 
- $X_t$ is a standard Brownian motion.
- $\gamma$ is the speed of reversion.
- $\bar{r}$ is the long-term mean level.
- $\sigma$ is the volatility.

Solution:

$$
r_t = r_0 e^{-\gamma t} + \bar{r}(1 - e^{-\gamma t}) + \sigma e^{-\gamma t} \int_0^t e^{\gamma s} d X_s
$$

Since $ \int_0^t e^{\gamma s} d X_s$ is normally distributed, $r_t$ is also normally distributed with mean and variance given by:

$$
E[r_t] = r_0 e^{-\gamma t} + \bar{r}(1 - e^{-\gamma t})
$$

$$
Var[r_t] = \frac{\sigma^2}{2\gamma}(1 - e^{-2\gamma t})
$$

### Simulating paths
#### Euler-Maruyama Method

$$
r_{t+\Delta t} = r_t + \gamma (\bar{r} - r_t) \Delta t + \sigma \sqrt{\Delta t} Z_t
$$

where $Z_t \sim N(0, 1)$

Euler-Maruyama is a straightforward numerical method for simulating SDEs. It approximates the continuous-time process by discrete steps, making it easy to implement and understand. However, it may require small time steps for accuracy. If $\Delta t$ is too large, the approximation may deviate significantly from the true process.

#### Transition Density Method

The transition density method leverages the known distribution of the Vasicek process (indicated by the mean and variance formulas above) to directly sample from the distribution at each time step. This method can be more efficient and accurate than Euler-Maruyama, which relies on discretizing an SDE.

$$
r_{t + \Delta t} = \bar{r} + (r_t - \bar{r}) e^{-\gamma \Delta t} + \sqrt{\frac{\sigma^2}{2\gamma}(1 - e^{-2\gamma \Delta t})} Z_t
$$


## Cox-Ingersoll-Ross (CIR) Model

SDE:

$$
d r_t = \gamma (\bar{r} - r_t) d t + \sqrt{\alpha r_t} d X_t
$$

The transition density of the CIR process is given by a non-central chi-squared distribution:

$$
f(r_{t+\Delta t} | r_t) = c_s \chi^2(c_s r_{t+\Delta t}; \nu, \lambda_{t+\Delta t})
$$

where

- $c_s = \frac{4\gamma}{\alpha(1 - e^{-\gamma \Delta t})}$
- $\nu = \frac{4\gamma \bar{r}}{\alpha}$
- $\lambda_{t+\Delta t} = c_s r_t e^{-\gamma \Delta t}$

which is much more complex than the Vasicek model (normally distributed). 

###
#### Euler-Maruyama Method

$$
r_{t+\Delta t} = r_t + \gamma (\bar{r} - r_t) \Delta t + \sqrt{\alpha r_t} \sqrt{\Delta t} Z_t
$$

However, this method can lead to negative interest rates if $r_t$ is small and the random shock is large enough.

#### Transition Density Method


$$
r_{t + \Delta t} = c_s \chi^2(c_s r_{t+\Delta t}; \nu, \lambda_{t+\Delta t})
$$

$$
\lambda_{t+\Delta t} = c_s r_t e^{-\gamma \Delta t}
$$

$$
\nu = \frac{4\gamma \bar{r}}{\alpha}
$$

