

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


### Calibration

#### Real World Calibration

We assume we have observed short rate data at discrete times: $r_0, r_{\Delta}, r_{2\Delta}, \ldots, r_{N\Delta}$, perhaps overnight rates observed over five years. 

##### Maximum Likelihood Estimation (MLE)

In the previous section, we derived that the transition density of the Vasicek process is normally distributed. Rewriting the mean and variance:



$$
E[r_{i \Delta} | r_{(i-1) \Delta}] = \alpha^* + \beta^* r_{(i-1) \Delta}
$$

$$
Var[r_{i \Delta} | r_{(i-1) \Delta}] = \sigma^{*2}, \quad i = 1, 2, \ldots, N
$$

where

$$
\alpha^* = \bar{r}(1 - e^{-\gamma \Delta})
$$

$$
\beta^* = e^{-\gamma \Delta}
$$

$$
\sigma^{*2} = \frac{\sigma^2}{2\gamma}(1 - e^{-2\gamma \Delta})
$$

The density function is:

$$
f(r_{i \Delta} | r_{(i-1) \Delta}) = \frac{1}{\sqrt{2\pi \sigma^{*2}}} \exp\left(-\frac{(r_{i \Delta} - (\alpha^* + \beta^* r_{(i-1) \Delta}))^2}{2\sigma^{*2}}\right)
$$

The joint likelihood function for the observed data is:

$$
L = f_0(r_0 | \alpha^*, \beta^*, \sigma^{*2}) \prod_{i=1}^{N} f(r_{i \Delta} | r_{(i-1) \Delta})\\
= f_0(r_0|\alpha^*, \beta^*, \sigma^{*2}) \exp\left(-\frac{1}{2\sigma^{*2}} \sum_{i=1}^{N} (r_{i \Delta} - (\alpha^* + \beta^* r_{(i-1) \Delta}))^2\right) \frac{1}{\sqrt{(2\pi)^N} \sigma^{*N}}
$$

If sample size $N$ is large, the influence of the initial density $f_0$ is negligible. 

Else, we can assume $r_0$ is normally distribution with mean $\frac{\alpha^*}{1 - \beta^*}$ and variance $\frac{\sigma^{*2}}{1 - \beta^{*2}}$. The pdf of $r_0$ is:

$$
f_0(r_0 | \alpha^*, \beta^*, \sigma^{*2}) = \frac{1}{\sqrt{2\pi \frac{\sigma^{*2}}{1 - \beta^{*2}}}} \exp\left(-\frac{(r_0 - \frac{\alpha^*}{1 - \beta^*})^2}{2 \frac{\sigma^{*2}}{1 - \beta^{*2}}}\right)
$$

Taking the logarithm of the likelihood function, we get the log-likelihood function:

$$
\ln L = - \frac{1}{2} ln(2 \pi) - \frac{1}{2} ln(2 \sigma^{*2} / (1 - \beta^{*2})) - \frac{(r_0 - [\alpha^*/(1 - \beta^*)])}{2 \sigma^{*2} / (1 - \beta^{*2})} - \frac{N}{2} ln(2 \pi) - N ln(\sigma^*) - \frac{1}{2\sigma^{*2}} \sum_{i=1}^{N} (r_{i \Delta} - (\alpha^* + \beta^* r_{(i-1) \Delta}))^2
$$


#### Risk Neutral Calibration

Under risk neutral calibration, we have a set of zero-coupon bond prices observed in the market for different maturities. Different from the real-world calibration, we need to adjust the parameters to their risk-neutral equivalents to price these bonds correctly rather than fitting historical short rate data.

With Vasicek model, we can derive the closed-form solution for zero-coupon bond prices:

$$
P(r, t, T) = e^{A(t, T) - B(t, T) r_t}
$$

where

$$
B(t, T) = \frac{1 - e^{-\gamma^* (T - t)}}{\gamma^*}
$$

$$
A(t, T) = (B(t, T) - (T - t))(\bar{r}^* - \frac{\sigma^2}{2\gamma^{*2}}) - \frac{\sigma^2}{4\gamma^*} B(t, T)^2
$$

The parameters $\gamma^*$, $\bar{r}^*$, and $\sigma$ are the risk-neutral parameters with is different from the real-world parameters $\gamma$, $\bar{r}$, and $\sigma$. 

To calibrate the model, we minimize the sum of squared differences between the market prices of zero-coupon bonds and the model prices:

$$
\min_{\gamma^*, \bar{r}^*, \sigma} \sum_{i=1}^{N} (P^{mkt}(0, T_i) - P(0, T_i))^2
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

### Simulating paths
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


### Calibration

#### The generalized method of moments (GMM)

To use GMM for the CKLS short rate model given by:

$$
d r_t = \gamma (\bar{r} - r_t) d t + (\alpha)^{1/2} r_t^{\tau} d X_t
$$

CIR model is a special case with $\tau = 0.5$.

Define the moment conditions based on the conditional mean and variance of the process:

$$
f_1 = \frac{1}{N} \sum_{i=1}^{N} (r_{i \Delta} - \alpha_1 - \beta_1 r_{(i-1) \Delta})
$$

$$
f_2 = \frac{1}{N} \sum_{i=1}^{N} ((r_{i \Delta} - \alpha_1 - \beta_1 r_{(i-1) \Delta})^2 - r_{(i-1) \Delta} \sigma^2)
$$

$$
f_3 = \frac{1}{N} \sum_{i=1}^{N} (r_{i \Delta} - \alpha_1 - \beta_1 r_{(i-1) \Delta}) r_{(i-1)\Delta}
$$

Choose $\alpha_1$, $\beta_1$, and $\sigma^2$ to minimize:

$$
f_1^2 + f_2^2 + f_3^2
$$

the estimates for the parameters are:

$$
\gamma = \frac{1 - \beta_1}{\Delta}
$$

$$
\bar{r} = \frac{\alpha_1}{1 - \beta_1}
$$

$$
\alpha = \frac{\sigma^2}{\Delta}
$$

#### Maximum Likelihood Estimation (MLE)

$$
L = f_0(r_0 | \gamma, \bar{r}, \alpha) \prod_{i=1}^{N} c_{\Delta} \chi^2(c_{\Delta} r_{i \Delta}; \nu, \lambda_{i \Delta})
$$

where
- $c_{\Delta} = \frac{4\gamma}{\alpha(1 - e^{-\gamma \Delta})}$
- $\nu = \frac{4\gamma \bar{r}}{\alpha}$
- $\lambda_{i \Delta} = c_{\Delta} r_{(i-1) \Delta} e^{-\gamma \Delta}$
- we can ignore $f_0$.

### Risk Neutral Calibration

With CIR model, we can derive the closed-form solution for zero-coupon bond prices:

$$
P(r, t, T) = e^{A(t, T) - B(t, T) r_t}
$$

where

$$
B(t, T) = \frac{2(e^{\eta (T - t)} - 1)}{(\gamma^* + \eta)(e^{\eta (T - t)} - 1) + 2\eta}
$$

$$
A(t, T) = \frac{2\gamma^* \bar{r}^*}{\alpha} \ln\left(\frac{2\eta e^{(\gamma^* + \eta)(T - t)/2}}{(\gamma^* + \eta)(e^{\eta (T - t)} - 1) + 2\eta}\right)
$$

$$
\eta = \sqrt{\gamma^{*2} + 2\alpha}
$$

minimizing the sum of squared differences between the market prices of zero-coupon bonds and the model prices given $\alpha$:

$$
\min_{\gamma^*, \bar{r}^*} \sum_{i=1}^{N} (P^{mkt}(0, T_i) - P(0, T_i))^2
$$

Note that $\alpha$ is estimated from the real-world data and kept constant during the risk-neutral calibration instead of estimated together with $\gamma^*$ and $\bar{r}^*$. This is because the volatility parameter $\alpha$ is assumed to be the same under both the real-world and risk-neutral measures, while the drift parameters $\gamma$ and $\bar{r}$ can differ between the two measures due to risk premia. Also, the influence of $\alpha$ on bond prices is typically less significant compared to the drift parameters, making the calibration unstable if all three parameters are estimated simultaneously.

$\gamma^* \times \bar{r}^*$ should be greater than $\frac{\alpha}{2}$ to ensure the Feller condition is satisfied and the interest rates remain positive.




## Two Factor Vasicek Model

The model is a short rate process as the sum of two factors, a long-term factor and a short-term factor:

$$
r_t = \phi_{1,t} + \phi_{2,t}
$$

where each factor follows its own Vasicek process:

$$
d \phi_{i, t} = \gamma_i (\bar{\phi}_i - \phi_{i,t}) d t + \sigma_i d X_{i,t}, \quad i = 1, 2
$$

$$
d X_{1,t} d X_{2,t} = \rho d t
$$

The solution of the SDEs for each factor is:

$$
\phi_{i, t + s} = \\phi_{i, t} e^{-\gamma_i s} + \bar{\phi}_i (1 - e^{-\gamma_i s}) + \sigma_i e^{-\gamma_i s} \int_t^{t+s} e^{\gamma_i u} d X_{i,u}, \quad i = 1, 2
$$

The mean and variance of each factor are:

$$
E[\phi_{i, t + s} | \phi_{i, t}] = \phi_{i, t} e^{-\gamma_i s} + \bar{\phi}_i (1 - e^{-\gamma_i s}), \quad i = 1, 2
$$

$$
Var[\phi_{i, t + s} | \phi_{i, t}] = \frac{\sigma_i^2}{2\gamma_i}(1 - e^{-2\gamma_i s}), \quad i = 1, 2
$$

Using Ito's isometry, the covariance between the two factors is:

$$
Cov[\phi_{1, t + s}, \phi_{2, t + s} | \phi_{1, t}, \phi_{2, t}] = \frac{\rho \sigma_1 \sigma_2}{\gamma_1 + \gamma_2}(1 - e^{-(\gamma_1 + \gamma_2)s})
$$

The conditional correlation:

$$
\rho(s) = Corr[\phi_{1, t + s}, \phi_{2, t + s} | \phi_{1, t}, \phi_{2, t}] = \frac{Cov[\phi_{1, t + s}, \phi_{2, t + s} | \phi_{1, t}, \phi_{2, t}]}{\sqrt{Var[\phi_{1, t + s} | \phi_{1, t}] Var[\phi_{2, t + s} | \phi_{2, t}]}} \\
= \rho (\frac{4 \gamma_1 \gamma_2 (1 - e^{-(\gamma_1 + \gamma_2)s})^2}{(\gamma_1 + \gamma_2)^2 (1 - e^{-2\gamma_1 s})(1 - e^{-2\gamma_2 s})})^{1/2}
$$

### Simulating paths - Transition Density Method

Steps:
1. Choose initial short term rate $r_0$ and long term rate $r_0(\tau)$ for a given maturity $\tau$.
2. Decompose for $\phi_{1,0}$ and $\phi_{2,0}$:
    
    $$
    \phi_{1,0} + \phi_{2,0} = r_0
    $$

    $$
    A(\tau) - B_1(\tau) \phi_{1,0} - B_2(\tau) \phi_{2,0} = - r_0(\tau) \tau
    $$

    $$
    \phi_{1,0} = \frac{r_0(\tau) \tau + A(\tau) - B_2(\tau) r_0}{B_1(\tau) - B_2(\tau)}
    $$

    $$
    \phi_{2,0} = r_0 - \phi_{1,0}
    $$

3. Simulate standard bivariate normal random variables $(Z_{1,t}, Z_{2,t})$ with correlation $\rho(s)$.
4. Update each factor using the transition density:
    
    $$
    \phi_{i, t+s} = \phi_{i, t} e^{-\gamma_i s} + \bar{\phi}_i (1 - e^{-\gamma_i s}) + \sqrt{\frac{\sigma_i^2}{2\gamma_i}(1 - e^{-2\gamma_i s})} Z_{i,t}, \quad i = 1, 2
    $$

    $$
    r_{t+s} = \phi_{1, t+s} + \phi_{2, t+s}
    $$
