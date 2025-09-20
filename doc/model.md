

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
- $r_t$ is the short rate at time $t$. 

Solution:

$$
r_{t+s} = r_t e^{-\gamma s} + \bar{r}(1 - e^{-\gamma s}) + \sigma e^{-\gamma s} \int_0^s e^{\gamma u} d X_u
$$

Since $ \int_0^s e^{\gamma u} d X_u$ is normally distributed, $r_t$ is also normally distributed with mean and variance given by:

$$
E[r_t, s] = \bar{r} + (r_t - \bar{r}) e^{-\gamma s}
$$


$$
\sigma^2 = \frac{\sigma^2}{2\gamma}(1 - e^{-2\gamma s})
$$

### Simulating paths
#### Euler-Maruyama Method

To simulate, we discretize the time into small intervals of size $\Delta t$ and update the short rate using the Euler-Maruyama approximation:

$$
r_{t+\Delta t} = r_t + \gamma (\bar{r} - r_t) \Delta t + \sigma \sqrt{\Delta t} Z_t
$$

where $Z_t \sim N(0, 1)$ and the initial short rate $r_0$ is normally overnight rate from the market.

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
\ln L = - \frac{1}{2} ln(2 \pi) - \frac{1}{2} ln(2 \sigma^{*2} / (1 - \beta^{*2})) - \frac{(r_0 - [\alpha^*/(1 - \beta^*)])^2}{2 \sigma^{*2} / (1 - \beta^{*2})} \\
- \frac{N}{2} ln(2 \pi) - N ln(\sigma^*) - \frac{1}{2\sigma^{*2}} \sum_{i=1}^{N} (r_{i \Delta} - (\alpha^* + \beta^* r_{(i-1) \Delta}))^2
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
r_{t + \Delta t} \sim \frac{1}{c_{\Delta t}} \chi^2(c_{\Delta t} r_{t+\Delta t}; \nu, \lambda_{t+\Delta t})
$$

$$
\lambda_{t+\Delta t} = c_{\Delta t} r_t e^{-\gamma \Delta t}
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

$\gamma^* \times \bar{r}^*$ should be greater than $\frac{\alpha}{2}$ to ensure the Feller condition is satisfied and the interest rates remain positive. If not, the model may not be suitable for the market data.




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

### Calibration

The real-world calibration of the two-factor model is more complex. Therefore, we only carlibrate the model under the risk-neutral measure using market prices of zero-coupon bonds.

Data needed:
- Volatility of short-term rate $\sigma_{\tau}$, which can be estimated from historical short rate data.
- Volatility of long-term rate $\sigma_{LT}$, which can be estimated from historical long-term rate data.
- Correlation $\rho$ between short-term and long-term rates, which can be estimated from historical data of both rates.
- Market prices of zero-coupon bonds for different maturities.


The closed-form solution for zero-coupon bond prices is:

$$
Z(\phi_{1,t}, \phi_{2,t}, t, T) = e^{A(t, T) - \phi_{1,t} B_1(t, T) - \phi_{2,t} B_2(t, T)}
$$

where

$$
B_i(t, T) = \int_0^{T-t} e^{-\gamma_i s} ds, i = 1, 2, 3
$$

$$
\begin{aligned}
A(t;T) &= 
\vec{\phi}_1 \big( B_1(t;T) - (T-t) \big) 
- \frac{\sigma_1^2}{2} \left[ 
\frac{1}{\gamma_1^2} \big( B_1(t;T) - (T-t) \big) 
+ \frac{1}{2\gamma_1} B_1(t;T)^2
\right] \\
&\quad + \vec{\phi}_2 \big( B_2(t;T) - (T-t) \big) 
- \frac{\sigma_2^2}{2} \left[ 
\frac{1}{\gamma_2^2} \big( B_2(t;T) - (T-t) \big) 
+ \frac{1}{2\gamma_2} B_2(t;T)^2
\right] \\
&\quad + \frac{\rho \sigma_1 \sigma_2}{\gamma_1 \gamma_2} 
\big[ B_1(t;T) + B_2(t;T) - B_3(t;T) - (T-t) \big]
\end{aligned}
$$

where

$$
\gamma_3 = \gamma_1 + \gamma_2
$$

We minimize the sum of squared differences between the market prices of zero-coupon bonds and the model prices:

$$
J(\gamma^*, \bar{r}^*) = \sum_{i=1}^{N} (P^{mkt}(0, T_i) - P(0, T_i))^2
$$



- **When \(\gamma_1 = 0, \ \gamma_2 \neq 0\)**  

$$
\begin{aligned}
A(t;T) &= 
\vec{\phi}_2 \big( B_2(t;T) - (T-t) \big) 
+ \sigma_1^2 \frac{(T-t)^3}{6} \\
&\quad - \frac{\sigma_2^2}{2} \left[ 
\frac{1}{\gamma_2^2} \big( B_2(t;T) - (T-t) \big) 
+ \frac{1}{2\gamma_2} B_2(t;T)^2
\right] \\
&\quad + \rho \sigma_1 \sigma_2 
\left[
\frac{(T-t)^2}{2} + \frac{\exp(-\gamma_2 (T-t))(\gamma_2 (T-t) + 1) - 1}{\gamma_2^2}
\right]
\end{aligned}
$$


- **When \(\gamma_1 \neq 0, \ \gamma_2 = 0\)**  

$$
\begin{aligned}
A(t;T) &= 
\vec{\phi}_1 \big( B_1(t;T) - (T-t) \big) 
- \frac{\sigma_1^2}{2} \left[ 
\frac{1}{\gamma_1^2} \big( B_1(t;T) - (T-t) \big) 
+ \frac{1}{2\gamma_1} B_1(t;T)^2
\right] \\
&\quad + \sigma_2^2 \frac{(T-t)^3}{6} \\
&\quad + \rho \sigma_1 \sigma_2 
\left[
\frac{(T-t)^2}{2} + \frac{\exp(-\gamma_1 (T-t))(\gamma_1 (T-t) + 1) - 1}{\gamma_1^2}
\right]
\end{aligned}
$$

- **When \(\gamma_1 = 0, \ \gamma_2 = 0\)**

$$
A(t; T) = (\sigma_1^2 + \sigma_2^2 + 2 \rho \sigma_1 \sigma_2) \frac{(T-t)^3}{6}
$$


In each step, we fix $\sigma_1$, $\sigma_2$, and $\rho$ to minimize $J(\gamma^*, \bar{r}^*)$. Then we fix $\gamma_1^*$, $\gamma_2^*$, and $\bar{r}^*$ and update $\sigma_1$, $\sigma_2$, and $\rho$ based on the functions:

$$
\sigma_{\tau} = \sqrt{\sigma_1^2 + \sigma_2^2 + 2 \rho \sigma_1 \sigma_2}
$$


$$
\sigma(\tau) = \sqrt{ 
\sigma_1^2 \left(\frac{B_1(\tau)}{\tau}\right)^2
+ \sigma_2^2 \left(\frac{B_2(\tau)}{\tau}\right)^2
+ 2 \left(\frac{B_1(\tau)}{\tau}\right)\left(\frac{B_2(\tau)}{\tau}\right)\sigma_1\sigma_2\rho
}
$$

$$
\rho(0,\tau) = \frac{
\sigma_1^2 \frac{B_1(\tau)}{\tau}
+ \sigma_2^2 \frac{B_2(\tau)}{\tau}
+ \left(\frac{B_1(\tau)}{\tau} + \frac{B_2(\tau)}{\tau}\right)\sigma_1\sigma_2\rho
}{
\sigma_\tau \, \sigma(\tau)
}
$$


Use 'nls' in R to estimate and calculate the p value of each parameter to check if it is significant. If not, this model may not be suitable for the market data.

# No-Arbitrage Term Structure Models

The equilibrium models are not exact fits, they just minimize the total errors between model prices and market prices. To eliminate the mismatch, Hull-White extend Vasicek models by adding a time-dependent function to the drift term.

In theory, no-arbitrage models reproduce the yield curve and prices for a set of interest rate derivatives such as caps, floors or swaptions. In the calibration process, we will see that they do not exactly fit the market prices, but the errors are very small.


## Hull-White One Factor Model

An extension of the Vasicek model with a time-dependent function $\theta(t)$ added to the drift term:

$$
d r_t = [\theta(t) - \gamma^* r_t] d t + \sigma d X_t
$$

The solution of the SDE is:

$$
r_{t+s} = r_t \exp(-\gamma^* s) 
+ \exp(-\gamma^*(t+s)) \int_t^{t+s} \theta_u \exp(\gamma^* u)\,du 
+ \sigma \exp(-\gamma^* s) \int_0^s \exp(\gamma^* u)\, dX_u
$$

$$
E[r_{t+s} \mid r_t] = r_t \exp(-\gamma^* s) 
+ \exp(-\gamma^*(t+s)) \int_t^{t+s} \theta_u \exp(\gamma^* u)\,du\\
= r_t \exp(-\gamma^* s) \\
+ f(0, s+t) - f(0,t)\exp(-\gamma^* s) \\
+ \frac{\sigma^2}{2(\gamma^*)^2} \Big[ 
1 - \exp(-\gamma^* s) 
+ \exp(-2\gamma^*(t+s)) 
- \exp(-\gamma^*(2t+s)) 
\Big]
$$

$$
Var[r_{t+s} \mid r_t] 
= \big[\sigma \exp(-\gamma^* s)\big]^2 \int_0^s \exp(2\gamma^* u)\,du
= \frac{\sigma^2}{2\gamma^*}\left(1 - e^{-2\gamma^* s}\right)
$$

Here $f(0,t)$ is the instantaneous forward rate at time 0 for maturity $t$ derived from the initial yield curve.

The Zero-coupon bond price is given by:

$$
Z(r_t, t; T) = e^{A(t;T) - B(t;T) r_t}
$$

$$
B(t;T) = \frac{1 - e^{-\gamma^* (T - t)}}{\gamma^*}
$$

$$
A(t;T) = - \int_t^T B(t, T-u) \theta(u)\,du \\
+ \frac{\sigma^2}{2 \gamma^{*2}} (T + \frac{1-e^{-2 \gamma^* (T-t)}}{2 \gamma^*} - 2 B(t;T)) \\ 
$$

By taking partial derivative of $Z(r_0, 0; T)$ with respect to $T$, we can derive $\theta(t)$ as a function of $f(0,t)$:

$$
f(0,t) = -\frac{\partial \log Z(r_0, 0; t)}{\partial t}
$$

which is the instantaneous forward rate at time 0 for maturity $t$ derived from the initial yield curve.

$$
\theta(t) = \frac{\partial f(0,t)}{\partial t} + \gamma^* f(0,t) + \frac{\sigma^2}{2\gamma^*}(1 - e^{-2\gamma^* t})
$$

### Simulating paths

#### Transition Density Method

As the Hull-White model is a Gaussian model (Ornstein-Uhlenbeck process), the short rate is normally distributed.
The expected value and variance of $r_{t+s}$ are given above. Therefore, we can simulate the short rate directly:

$$
r_{t + \Delta t} = E[r_{t + \Delta t} | r_t] + \sqrt{Var[r_{t + \Delta t} | r_t]} Z_t
$$

#### Euler-Maruyama Method

Similarly, we can use Euler-Maruyama method to simulate the short rate:

$$
r_{t+\Delta t} = r_t + [\theta(t) - \gamma^* r_t] \Delta t + \sigma \sqrt{\Delta t} Z_t
$$

### Calibration

We should firstly interpolate the yield curve with cubic splines, higher degree polynomials or Nelson-Siegel method to get a smooth curve. $f(0,t)$ can be derived from the fitted yield curve.

(extra step: compare the interpolation methods. If the forward rates are oscillating, try another method.)

Steps:
1. obtain spot rates $r(0,t)$ from the yield curve for different maturities $t$.
2. Fit a curve, $\hat{r}(0, t)$
3. Calculate $f(0,t)$ from $\hat{r}(0,t)$:

$$
\hat{f}(0,t) = \frac{\partial t \hat{r}(0,t)}{\partial t}
$$

Use the data table 19.4 which gives swap rate, discount factors and cap prices to calibrate the model.

Steps in calibration:
1. Estimate forward rates $f(0,t)$ from the yield curve.
2. Calibrate implied volatilities of caps.
3. Choose initial guess for $\gamma^*$ and $\sigma$ and calculate $\theta(t)$.
4. Calibrate $\gamma^*$ and $\sigma$ by minimizing the sum of squared differences between market prices of caps and model prices of caps (or implied volatilities of caps)
5. Repeat steps 2 and 3 until convergence and obtain $\theta(t)$ using $f(0,t)$, $\gamma^*$ and $\sigma$.

---

The price of a caplet is given by the Black formula:

$$
V(r_0,0) = M \times \Big( K Z(r_0,0;T-\Delta)N(-d_2) - Z(r_0,0;T)N(-d_1) \Big)
$$

$$
d_1 = \frac{1}{S_Z(T-\Delta;T)} \log\!\left( \frac{Z(r_0,0;T)}{K Z(r_0,0;T-\Delta)} \right) 
+ \frac{S_Z(T-\Delta;T)}{2}
$$

$$
d_2 = d_1 - S_Z(T-\Delta;T)
$$

$B(t;T)$ is given by:

$$
B(t;T) = \frac{1 - e^{-\gamma^* (T - t)}}{\gamma^*}
$$

$Z(r_0,0;T)$ is the zero-coupon bond price at time 0 for maturity $T$


To price caps, we need to calculate the volatility of the zero-coupon bond price:

$$
S_Z(T_o;T_B)^2 = B(T_o;T_B)^2 \cdot \frac{\sigma^2}{2\gamma}\left(1 - e^{-2\gamma T_o}\right)
$$
where $T_o$ is the option maturity or reset date, $T_B$ is the bond maturity which is the caplet payment date.

$$
A(t;T) = \log\!\left(\frac{Z(r_0,0;T)}{Z(r_0,0;t)}\right) 
+ B(t;T) f(0,t) - \frac{\sigma^2}{4\gamma} B(t;T)^2 \left(1 - e^{-2\gamma t}\right)
$$

$$
A(t;T) = \log\!\left(\frac{Z(r_0,0;T)}{Z(r_0,0;t)}\right) 
+ (T-t)f(0,t) - \frac{\sigma^2}{2}(T-t)^2 t
$$

$$
K = \frac{1}{1+r_K\Delta}
$$

Notice that $\theta(t)$ is not needed in the caplet pricing formula. That's because $\theta(t)$ is determined by the initial yield curve and does not affect the dynamics of the short rate under the risk-neutral measure. In derivatives pricing, we are primarily concerned with the volatility and mean-reversion characteristics of the short rate, which are captured by $\sigma$ and $\gamma^*$. Therefore, only these two parameters need to be calibrated to market prices of interest rate derivatives.


## Two Factor Hull-White Model

An extension of the two-factor Vasicek model with time-dependent functions $\theta_1(t)$ and $\theta_2(t)$ added to the drift terms:

$$
r_t = \phi_{1,t} + \phi_{2,t}
$$

where each factor follows its own Hull-White process:

$$
d \phi_{1, t} = (\theta_t - \gamma_1^* \phi_{1,t}) d t + \sigma_1 d X_{1,t}
$$

$$
d \phi_{2, t} = - \gamma_2^* \phi_{2,t} d t + \sigma_2 d X_{2,t}
$$

$$
d X_{1,t} d X_{2,t} = \rho d t
$$


The solution of the SDEs for each factor is:

$$
\phi_{1,t+s} = \phi_{1,t}\exp(-\gamma_1^* s) 
+ \exp(-\gamma_1^*(t+s)) \int_t^{t+s} \theta_u \exp(\gamma_1^* u)\,du 
+ \sigma_1 \exp(-\gamma_1^* s) \int_0^s \exp(\gamma_1^* u)\, dX_{1,u}
$$

$$
\phi_{2,t+s} = \phi_{2,t}\exp(-\gamma_2^* s) 
+ \sigma_2 \exp(-\gamma_2 s) \int_0^s \exp(\gamma_2 u)\, dX_{2,u}
$$


$$
E[\phi_{1,t+s} \mid \phi_{1,t}] 
= \phi_{1,t}\exp(-\gamma_1^* s) 
+ \exp(-\gamma_1^*(t+s)) \int_t^{t+s} \theta_u \exp(\gamma_1^* u)\,du
$$

$$
E[\phi_{2,t+s} \mid \phi_{2,t}] 
= \phi_{2,t}\exp(-\gamma_2 s)
$$

$$
Var[\phi_{1,t+s} \mid \phi_{1,t}] 
= \sigma_1^2 \frac{1 - \exp(-2\gamma_1^* s)}{2\gamma_1^*}
$$

$$
Var[\phi_{2,t+s} \mid \phi_{2,t}] 
= \sigma_2^2 \frac{1 - \exp(-2\gamma_2 s)}{2\gamma_2}
$$

$$
Cov[\phi_{1,t+s}, \phi_{2,t+s} \mid \phi_{1,t}, \phi_{2,t}] 
= \frac{\rho\sigma_1\sigma_2}{\gamma_1+\gamma_2} 
\Big( 1 - \exp(-( \gamma_1 + \gamma_2 )s) \Big)
$$

### Calibration

**Forward rate volatility term (Two-Factor Hull–White)**
$$
S_Z(T_o;T_B)^2
= B_1(T_o;T_B)^2 \,\frac{\sigma_1^2}{2\gamma_1^*}\big(1-e^{-2\gamma_1^* T_o}\big) \\
+ B_2(T_o;T_B)^2 \,\frac{\sigma_2^2}{2\gamma_2^*}\big(1-e^{-2\gamma_2^* T_o}\big) \\
+ B_1(T_o;T_B) B_2(T_o;T_B)\, \frac{1 - e^{-(\gamma_1^*+\gamma_2^*) T_o}}{\gamma_1^*+\gamma_2^*}\,\sigma_1\sigma_2\rho
$$

where
- $T_o$: reset date, when the forward rate for the period is fixed;  
- $T_B$: payment date, when the cash flow is settled;  
- $\Delta = T_B - T_o$ is the accrual period length;  
- $B_i(t;T) = \frac{1 - e^{-\gamma_i^*(T-t)}}{\gamma_i^*}$ is the standard factor in the two-factor Hull–White/Vasicek model.  

**Strike re-expression**
$$
K = \frac{1}{1 + r_K \Delta}
$$

**Caplet value (Hull–White Gaussian closed form)**
$$
V(r_0,0) = M \,\Big( K\,Z(r_0,0;T-\Delta)\,N(-d_2)\;-\; Z(r_0,0;T)\,N(-d_1)\Big)
$$

Given the market prices of caps, we need to bootstrap the caplet prices first. 

with
$$
d_1 = \frac{1}{S_Z(T-\Delta;T)}\,
\log\!\left(\frac{Z(r_0,0;T)}{K\,Z(r_0,0;T-\Delta)}\right)
+ \frac{S_Z(T-\Delta;T)}{2},
\qquad
d_2 = d_1 - S_Z(T-\Delta;T).
$$

However, the calibration is very sensitive to the initial guess of the parameters and may find a local minimum. Instead, we can use implied volatilities from caplet prices to calibrate the model. 

To do this, we need:
1. Bootstrap the caplet prices from the market cap prices.
2. Calculate the implied volatilities from the caplet prices using the Black formula.
3. Calibrate the model parameters to minimize the sum of squared differences between the model-implied volatilities and market-implied volatilities. The model-implied volatilities are given by:

$$
S_Z(T-\Delta;T)^2 =
B_1(T-\Delta;T)^2 \cdot \frac{\sigma_1^2}{2\gamma_1^*}\big(1 - e^{-2\gamma_1^*(T-\Delta)}\big) \\ 
+ B_2(T-\Delta;T)^2 \cdot \frac{\sigma_2^2}{2\gamma_2^*}\big(1 - e^{-2\gamma_2^*(T-\Delta)}\big) \\
+ B_1(T-\Delta;T) B_2(T-\Delta;T) \cdot 
\frac{1 - e^{-(\gamma_1^*+\gamma_2^*)(T-\Delta)}}{\gamma_1^*+\gamma_2^*}\,\sigma_1\sigma_2\rho.
$$

This links the caplet payoff variance $S_Z$ directly to the model parameters, making calibration more stable.

For the one-factor Hull–White model, the formula simplifies to:

$$
S_Z(T-\Delta;T)^2 
= B(0;\Delta)^2 \cdot \frac{\sigma^2}{2\gamma^*}\Big(1 - e^{-2\gamma^*(T-\Delta)}\Big),
$$

where  

$$
B(0;\Delta) = \frac{1 - e^{-\gamma^*\Delta}}{\gamma^*}.
$$


# Applications

The whole process of modeling interest rates can be summarized in the following steps:
1. Build models based on the above theories, which should include simulation and calibration functions. A typical structure is:
    - __init__ function: initialize the model with parameters.
    - simulate function: simulate paths of short rates.
    - calibrate function: calibrate the model parameters to market data.


Model validation is also important to ensure the model is reliable and robust. The validation process includes:

1. Simulate some parameters to validate the models. If the models cannot recover the parameters well, the model calibration may not be reliable.
    - simulate parameters
    - simulate paths
    - calibrate the model
    - calculate:
        - bias
        - SD
        - RMSE
2. Apply the models to real market data and reprice the caps to check the model performance.
3. Out-of-sample test (Depend on data): use the model to price other interest rate derivatives such as swaptions and compare with market prices. Or simulate paths then compare the mean or median of the simulated rates with the actual rates in the future.
4. Sensitivity analysis: change the initial parameters or length of data or interpolation methods to check the stability of the models.
5. Styly analysos: 
    - long term rate is higher than short term rate
    - volatility of short term rate is higher than long term rate
    - correlation between short term and long term rate is positive

 
