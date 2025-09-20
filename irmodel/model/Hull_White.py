"""
Hull-White model implementation in Python 

"""

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq, minimize
from scipy.interpolate import CubicSpline
@dataclass
class HullWhiteParams:
    Ts: np.ndarray     # Maturity times for thetas
    thetas: np.ndarray # Mean reversion level
    gamma: float      # Mean reversion speed
    sigma: float      # Volatility


class HullWhiteCapsPricer:
    
    """
    Hull-White model implementation for cap pricing.
    """

    def d1(self, vol, z_start, z_end, strike_price) -> float:
        """
        Calculate d1 for the Black formula.

        Args:
            vol (float): Volatility.
            z_start (float): Start time.
            z_end (float): End time.
            strike_price (float): Strike price.
        """
        return 1 / vol * np.log(z_end / (z_start * strike_price)) + 0.5 * vol

    def d2(self, vol, z_start, z_end, strike_price) -> float:
        """
        Calculate d2 for the Black formula.

        Args:
            vol (float): Volatility.
            z_start (float): Start time.
            z_end (float): End time.
            strike_price (float): Strike price.
        """
        return self.d1(vol, z_start, z_end, strike_price) - vol


    def caplet_price(self, 
                     vol: float,
                      strike_price: float,
                      notional: float,
                      z_start: float,
                      z_end: float) -> float:
        """
        Price a cap using the closed-formula under the Hull-White model.

        Args:
            vol (float): Volatility.
            strike_price (float): Strike price of the cap.
            notional (float): Notional amount of the cap.
            z_start (float): discount factor at the start time of the caplet.         
            z_end (float): discount factor at the end time of the caplet

        """

        return notional * (strike_price * z_start * norm.cdf(- self.d2(vol, z_start, z_end, strike_price)) - 
                          z_end * norm.cdf(- self.d1(vol, z_start, z_end, strike_price)))

    def cap_price(self,
                  strike_price: float,
                  vols: np.ndarray,
                  Ts: np.ndarray,
                  notional: float,
                  DFs: np.ndarray
                  ) -> float:
        """
        Price a cap using the closed-formula under the Hull-White model.

        Args:
            strike_price (float): Strike price of the cap.
            vols (np.ndarray): Volatilities for each caplet. shape = (len(Ts)-1,)
            Ts (np.ndarray): Maturity times of each caplet. The first element is the start time t which is not a maturity time.
            notional (float): Notional amount of the cap.
            DFs (np.ndarray): Discount factors for T. shape = (len(Ts),), The first element is the discount factor at time t.
        """

        if len(vols) >= len(Ts):
            raise ValueError("Length of Ts must be greater than length of vols.")

        price = 0.0
        for i in range(len(vols)):

            price += self.caplet_price(vols[i], strike_price, 
                                       notional, DFs[i], DFs[i + 1])

        return price


    def cal_implied_vol(self, 
                        market_price: float,
                        strike_price: float,
                        vols: np.ndarray,
                        Ts: np.ndarray,
                        notional: float,
                        DFs: np.ndarray,
                        initial_guess: float = 1,
                        ) -> dict:
        """
        calculate the implied volatility 

        Args:
            market_price (float): Market price of the cap.
            strike_price (float): Strike price of the cap.
            vols (np.ndarray): vols for caplets except the last one. shape = (len(Ts)-2,)
            Ts (np.ndarray): Maturity times of each caplet. The first element is the start time t which is not used.
            notional (float): Notional amount of the cap.
            DFs (np.ndarray): Discount factors for T. shape = (len(Ts),)
            initial_guess (float): Initial guess for volatility.
        """
        
        cap_price = self.cap_price(strike_price, vols, Ts[:-1], notional, DFs[:-1])
        caplet_price = market_price - cap_price

        def f(vol):
            return self.caplet_price(vol, strike_price, notional, 
                                      DFs[-2], DFs[-1]) - caplet_price
            
        result = brentq(f, 1e-6, 5.0, xtol=1e-6, maxiter=100)

        return {'implied_vol': result}

    def bootstrap_vols(self,
                       market_prices: np.ndarray,
                       strike_prices: np.ndarray,
                       Ts: np.ndarray,
                       notional: float,
                       DFs: np.ndarray,
                       initial_guess: float = 0.1
                       ) -> dict:
        """
        Bootstrap vols from market cap prices.

        Args:
            market_prices (np.ndarray): Market prices of the caps. shape = (n, )
            strike_prices (np.ndarray): Strike prices of the caps. shape = (n, )
            Ts (np.ndarray): Maturity times of each caplet. The first element is the start time t which is not used. shape = (n+1, )
            notional (float): Notional amount of the caps.
            DFs (np.ndarray): Discount factors for T. shape = (n+1, )
            initial_guess (float): Initial guess for volatility.
        """

        vols = np.zeros(len(market_prices))
        residual = 0

        for i in range(len(market_prices)):

            res = self.cal_implied_vol(market_prices[i], strike_prices[i], 
                                       vols[:i], Ts[:i+2], notional, DFs[:i+2], initial_guess)
            vols[i] = res['implied_vol']
            residual += (self.cap_price(strike_prices[i], vols[:i+1],
                                        Ts[:i+2], notional, DFs[:i+2]) - market_prices[i]) ** 2

        return {'success': True, 'vols': vols, 'error': residual,
                'rmse': np.sqrt(residual / len(market_prices))}


class f:
    
    def __f_and_df_cubic(self, 
                        rates: np.ndarray,
                        Ts: np.ndarray,
                        target_Ts: np.ndarray
                         ) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate forward rates and the differential of forward rates using polynomial interpolation.
        Args:
            rates (np.ndarray): Zero-coupon rates. shape = (n, )
            Ts (np.ndarray): Maturity times. shape = (n, )
            target_Ts (np.ndarray): Target maturity times to calculate forward rates. shape = (m, )
        Returns:
            tuple[np.ndarray, np.ndarray]: Forward rates and their differentials at target_Ts. shape = (m, ), (m, )
        """
        spl = CubicSpline(Ts, rates, bc_type='natural')
        r1 = spl.derivative(1)(target_Ts)
        r2 = spl.derivative(2)(target_Ts)
        r = spl(target_Ts)
        f = r + target_Ts * r1
        df = 2 * r1 + target_Ts * r2
        return f, df
        

    def __f_and_df_poly(self,
                        rates: np.ndarray,
                        Ts: np.ndarray,
                        target_Ts: np.ndarray,
                        n: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate forward rates and the differential of forward rates using polynomial interpolation.
        Args:
            rates (np.ndarray): Zero-coupon rates. shape = (n, )
            Ts (np.ndarray): Maturity times. shape = (n, )
            n (int): Degree of the polynomial.
        """
        p = np.polyfit(Ts, rates, n)
        r1 = np.polyder(p, 1)
        r2 = np.polyder(p, 2)
        r = np.polyval(p, target_Ts)
        f = r + target_Ts * np.polyval(r1, target_Ts)
        df = 2 * np.polyval(r1, target_Ts) + target_Ts * np.polyval(r2, target_Ts)
        return f, df


    def f_and_df(self, 
                 rates: np.ndarray,
                 Ts: np.ndarray,
                 target_Ts: np.ndarray,
                 n: int | None = None):
        """
        Calculate forward rates and the differential of forward rates.

        Args:
            rates (np.ndarray): Zero-coupon rates. shape = (n, )
            Ts (np.ndarray): Maturity times. shape = (n, )
            interpolation_method (str): Interpolation method. Options are 
                "poly" for polynomial interpolation (should fill n, e.g. n=3 for cubic), 
                "Nelson-Siegel" for Nelson-Siegel model,
        """
        if n == None:
            raise ValueError("Degree of polynomial n must be specified for polynomial interpolation.")
        elif n == 3:
            return self.__f_and_df_cubic(rates, Ts, target_Ts)
        else:
            return self.__f_and_df_poly(rates, Ts, target_Ts, n)

class HullWhiteFunc:
    
    @staticmethod
    def thetas(fs: np.ndarray, dfs: np.ndarray, gamma: float, sigma: float, Ts: np.ndarray) -> np.ndarray:
        """
        Calculate the mean reversion level theta(t) at each time step.
        """
        thetas = dfs + gamma * fs + sigma ** 2 / (2 * gamma) * (1 - np.exp(-2 * gamma * Ts))
        return thetas

    @staticmethod
    def B(gamma: float,
          T: float) -> float:
        """
        Calculate the function B(t, T) in the Hull-White model.
        """
        return (1 - np.exp(-gamma * T)) / gamma

    @staticmethod
    def A(gamma: float,
          sigma: float, 
          T: float,
          thetas: np.ndarray,
          theta_times: np.ndarray) -> float:
        """
        Calculate the function A(t, T) in the Hull-White model.
        """
        integral = 0.0
        for i in range(1, len(theta_times)):
            dt = min(T, theta_times[i]) - min(T, theta_times[i - 1])
            if dt > 0:
                integral += thetas[i - 1] * dt * HullWhiteFunc.B(gamma, T - theta_times[i - 1])

        return - integral + sigma ** 2 / (4 * gamma ** 2) * (T + (1 - np.exp(-2 * gamma * T)) / (2 * gamma) - 2 * HullWhiteFunc.B(gamma, T))

    @staticmethod
    def Z(gamma: float,
          sigma: float,
          T: float,
          r_0: float,
          thetas: np.ndarray,
          theta_times: np.ndarray) -> float:
        """
        Calculate the zero-coupon bond price Z(0, T) in the Hull-White model.
        """
        B = HullWhiteFunc.B(gamma, T)
        A = HullWhiteFunc.A(gamma, sigma, T, thetas, theta_times)
        return np.exp(A - B * r_0)




class HullWhiteCalibrator:

    def calibrate(self,
                  DFs: np.ndarray,
                  caps: np.ndarray,
                  cap_strikes: np.ndarray,
                  cap_maturities: np.ndarray,
                  notional: float = 100.0,
                  initial_guess: tuple[float, float] = (0.1, 0.01)
                  ) -> dict:
        """
        Calibrate Hull-White model parameters to market cap prices.
        Args:
            DFs (np.ndarray): Discount factors at each maturity time T and the start time t. shape = (n+1, )
            caps (np.ndarray): Market prices of the caps. shape = (n, )
            cap_strikes (np.ndarray): Strike prices of the caps. shape = (n, )
            cap_maturities (np.ndarray): Maturity times of the caps. shape = (n + 1, ). The first element is the start time t which is not used.
            notional (float): Notional amount of the caps.
        Returns:
            dict: Calibrated Hull-White gamma and sigma parameters, and optimization result.
        """
        
        # Calibrate implied vols from market cap prices
        pricer = HullWhiteCapsPricer()
        vol_res = pricer.bootstrap_vols(caps, cap_strikes, cap_maturities, notional, DFs) 

        if not vol_res['success']:
            raise ValueError("Vol bootstrapping failed.")

        # optimize gamma and sigma
        def objective(params):
            gamma, sigma = params
            if gamma <= 0 or sigma <= 0:
                return 1e10
            Bs = np.array([HullWhiteFunc.B(gamma, cap_maturities[i] - cap_maturities[i - 1]) \
                    for i in range(1, len(cap_maturities))])
            model_vols = Bs * sigma / np.sqrt(2 * gamma) * np.sqrt(1 - np.exp(-2 * gamma * cap_maturities[:-1]))
            return np.sum((model_vols - vol_res['vols']) ** 2)

        result = minimize(objective, initial_guess, bounds=[(1e-6, None), (1e-6, None)])
        
        return {"params": HullWhiteParams(
                    Ts=None,
                    thetas=None,
                    gamma=result.x[0],
                    sigma=result.x[1]),
                "success": result.success,
                "rmse": np.sqrt(result.fun / len(caps)),
                "rrmse": np.sqrt(result.fun) / len(caps) / np.mean(vol_res['vols']),
                }

class HullWhiteSimulator:
    def simulate(self, 
                 params: HullWhiteParams,
                 DFs: np.ndarray,
                 Ts: np.ndarray,
                 N: int,
                 dt: float = 1/252,
                 seed: int | None = None
                 ) -> np.ndarray:
        """
        Simulate short rate paths using the Hull-White model.
        Args:
            params (HullWhiteParams): Hull-White model parameters.
            DFs (np.ndarray): Discount factors at each maturity time T and the start time t. shape = (n, )
            Ts (np.ndarray): Maturity times. shape = (n, )
            n_scenarios (int): Number of scenarios to simulate.
            dt (float): Time step for simulation.
            seed (int | None): Random seed for reproducibility.
        Returns:
            np.ndarray: Simulated short rate paths. shape = (n_scenarios, n_time_steps)
        """
        if seed is not None:
            np.random.seed(seed)
        n_steps = int(Ts[-1] / dt) + 1
        short_rates = np.zeros((N, n_steps))

        gamma = params.gamma
        sigma = params.sigma
        thetas = params.thetas
        Ts = params.Ts
        for i in range(1, n_steps):
            dt_i = dt
            z = np.random.normal(0, 1, N)
            short_rates[:, i] = short_rates[:, i - 1] + (thetas[i - 1] - gamma * short_rates[:, i - 1]) * dt_i + sigma * np.sqrt(dt_i) * z

        return short_rates
        

