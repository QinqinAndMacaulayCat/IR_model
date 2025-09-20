"""
This module defines the Vasicek model for interatest rates, which is a type of one-factor short rate model.

Includes:
    - VasicekParams: A dataclass to hold the parameters of the Vasicek model.
    - Vasicek: A class that implements the Vasicek model, inheriting from BaseModel.
        - __init__: Initializes the Vasicek model with given parameters.
        - params: Property to access the model parameters.
        - simulate: Method to simulate the short rate path using the Vasicek model. use the Euler-Maruyama method or Transition Density method.
        - bond_price: Method to calculate the price of a zero-coupon bond under the Vasicek model.
        - calibrate: Method to calibrate the model parameters to market data. including Maximum Likelihood Estimation or Risk-Neutral Calibration.

References:
    doc-model.md
"""

import numpy as np
from scipy.optimize import least_squares, minimize
from dataclasses import dataclass, asdict

from .BaseModel import BaseModel


@dataclass
class VasicekParams:
    """
    Parameters for the Vasicek model.
    """
    gamma: float # Speed of reversion
    barR: float # Long-term mean
    sigma: float # Volatility

    def __post_init__(self):
        if self.sigma < 0:
            raise ValueError("Volatility must be positive.")


class VasicekSimulator:
    def simulate(self, 
                 params: VasicekParams,
                 r0: float, 
                 T: float, 
                 dt: float, 
                 N: int, 
                 method: str = "euler") -> np.ndarray:

        """
        Simulates the short rate path using the Vasicek model.
        Args:
            params (VasicekParams): Parameters of the Vasicek model.
            r0 (float): Initial short rate.
            T (float): Total time horizon.
            dt (float): Time step for simulation.
            N (int): Number of simulation paths.
            method (str): Method for simulation ("euler" or "transition").
        Returns:
            rates (np.ndarray): Simulated short rate path with shape (N+1,), the first element is r0.            
        """



        if method == "euler":
            return self.__simulate_euler(params, r0, T, dt, N) 
        elif method == "transition":
            return self.__simulate_transition(params, r0, T, dt, N)
        else:
            raise ValueError("Invalid method. Choose 'euler' or 'transition'.")

    def __simulate_euler(self, 
                        params: VasicekParams,
                        r0: float, T: float, dt: float, N: int) -> np.ndarray:
        """
        Simulates the short rate path using the Euler-Maruyama method.
        Args:
            params (VasicekParams): Parameters of the Vasicek model.
            r0 (float): Initial short rate.
            T (float): Total time horizon.
            dt (float): Time step for simulation.
            N (int): Number of simulation paths.
        Returns:
            rates (np.ndarray): Simulated short rate path with shape (int(T / dt) + 1, N), the first element is r0.
        """
        steps = int(T / dt)
        rates = np.zeros((steps + 1, N))
        rates[0] = r0
        gamma, barR, sigma = params.gamma, params.barR, params.sigma
        dW = np.random.normal(size = (steps, N)) * np.sqrt(dt) * sigma
        for i in range(1, steps + 1):
            dr = gamma * (barR - rates[i - 1]) * dt + dW[i - 1]
            rates[i] = rates[i - 1] + dr
        return rates

    def __simulate_transition(self, 
                              params: VasicekParams,
                              r0: float, 
                              T: float, 
                              dt: float, 
                              N: int) -> np.ndarray:
        """
        Simulates the short rate path using the Transition Density method.
        Args:
            params (VasicekParams): Parameters of the Vasicek model.
            r0 (float): Initial short rate.
            T (float): Total time horizon.
            dt (float): Time step for simulation.
            N (int): Number of simulation paths.
        Returns:
            rates (np.ndarray): Simulated short rate path with shape (steps + 1, N), the first element is r0.
        """
        steps = int(T / dt)
        rates = np.zeros((steps + 1, N))
        rates[0] = r0
        gamma, barR, sigma = params.gamma, params.barR, params.sigma
        z = np.random.normal(size = (steps, N))
        vol = sigma * np.sqrt((1 - np.exp(-2 * gamma * dt)) / (2 * gamma))


        for i in range(1, steps + 1):
            rates[i] = barR + (rates[i - 1] - barR) * np.exp(-gamma * dt) + vol * z[i - 1]

        return rates

class VasicekBondPricer:

    @ staticmethod
    def AB(gamma, barR, sigma, T) -> tuple[float, float]:
        """
        Computes the A and B functions for the Vasicek model bond pricing formula.
        """
        B = (1 - np.exp(- gamma * T)) / gamma    
        A = (B - T) * (barR - sigma**2 / (2 * gamma**2)) - (sigma**2 * B**2) / (4 * gamma)
        return (A, B)

    @staticmethod
    def risk_neutral_price(
                         r0: float, 
                         T: float, 
                         gamma: float, 
                         barR: float, 
                         sigma: float) -> float:    
        """
        Closed-form solution for zero-coupon bond price under the Vasicek model with face value 1.
        Args:
            r0 (float): Initial short rate.
            T (float): Time to maturity.
            gamma (float): Speed of reversion.
            barR (float): Long-term mean.
            sigma (float): Volatility.
        Returns:
            P (float): Price of the zero-coupon bond.
        """
        A, B = VasicekBondPricer.AB(gamma, barR, sigma, T)
        return np.exp(A - B * r0)




class VasicekCalibrator:

    def calibrate(self, 
                  method: str, 
                  *,
                  initial_params: VasicekParams,
                  r0: float | None = None,
                  zcb_prices: np.ndarray | None = None, 
                  maturities: np.ndarray | None= None, 
                  rates: np.ndarray | None = None, 
                  dt: float | None = None) -> dict:
        """
        Calibrates the model parameters to market data.
        Args:
            method (str): Calibration method ("mle" or "risk_neutral").
            initial_params (VasicekParams): Initial guess for the parameters.
            zcb_prices (np.ndarray): Market prices of zero-coupon bonds with face value 1 , r0 is not included.(required for "risk_neutral").
            maturities (np.ndarray): Corresponding maturities of the bonds (required for "risk_neutral").
            r0 (float): Initial short rate (required for "risk_neutral").
            rates (np.ndarray): Observed short rates (required for "mle").
            dt (float): Time step between observed rates (required for "mle").
        Returns:
            result (dict): Dictionary containing calibrated parameters and RMSE.
        """

        if method == "risk_neutral":
            if zcb_prices is None or maturities is None or r0 is None:
                raise ValueError("zcb_prices, maturities, and r0 are required for risk_neutral calibration.")
            return self.__calibrate_risk_neutral(initial_params, r0, zcb_prices, maturities)
        elif method == "mle":
            if rates is None or dt is None or r0 is None:
                raise ValueError("rates, dt, and r0 are required for mle calibration.")
            return self.__calibrate_mle(initial_params, r0, rates, dt)
        else:
            raise ValueError("Invalid method. Choose 'mle' or 'risk_neutral'.")

         
    def __calibrate_risk_neutral(self, 
                                 initial_params: VasicekParams,
                                 r0: float,
                                 zcb_prices: np.ndarray, 
                                 maturities: np.ndarray) -> dict:
        """
        Calibrates the model parameters to market data using Risk-Neutral Calibration.
        Args:
            initial_params (VasicekParams): Initial guess for the parameters.
            zcb_prices (np.ndarray): Market prices of zero-coupon bonds with face value 1.
            maturities (np.ndarray): Corresponding maturities of the bonds.
        Returns:
            result (dict): Dictionary containing calibrated parameters, and success flag.
        """

        def res(gamma, barR, sigma):
            model_prices = np.array([VasicekBondPricer.risk_neutral_price(r0, T, gamma, barR, sigma) for T in maturities])
            return model_prices - zcb_prices

        initial_guess = [initial_params.gamma, initial_params.barR, initial_params.sigma]
        bounds = ([1e-5, -np.inf, 1e-5], [np.inf, np.inf, np.inf])  # gamma > 0, sigma > 0
        result = least_squares(lambda x: res(x[0], x[1], x[2]), initial_guess, bounds=bounds)
        calibrated_params = VasicekParams(gamma=result.x[0], barR=result.x[1], sigma=result.x[2])
        return {"result": calibrated_params, 
                "success": result.success}

    def __likelihood(self, gamma: float, 
                     barR: float,
                     sigma: float,
                     r0: float,
                     short_rates: np.ndarray,
                     dt: float) -> float:
        """
        Computes the log-likelihood of the observed short rates given the model parameters.
        Args:
            gamma (float): Speed of reversion.
            barR (float): Long-term mean.
            sigma (float): Volatility.
            r0 (float): Initial short rate.
            short_rates (np.ndarray): Observed short rates (not including r0).
            dt (float): Time step between observed rates.
        Returns:
            log_likelihood (float): Log-likelihood of the observed short rates.
        """
        alpha = barR * (1 - np.exp(- gamma * dt))
        beta = np.exp(- gamma * dt)
        vol_sq = sigma**2 / (2 * gamma) * (1 - np.exp(- 2 * gamma * dt))
        a = - 0.5 * np.log(2 * np.pi) - 0.5 * np.log(2 * vol_sq / (1 - beta**2))
        b = - (r0 - alpha/(1 - beta))**2 / (2 * vol_sq / (1 - beta**2)) 
        c = - len(short_rates) / 2 * np.log(2 * np.pi) - len(short_rates) * 0.5 * np.log(vol_sq)
        prev = np.insert(short_rates[:-1], 0, r0)
        d = - np.sum((short_rates - alpha - beta * prev)**2)/ (2 * vol_sq)
        return a + b + c + d


    def __calibrate_mle(self, 
                        initial_params: VasicekParams, 
                        r0: float,
                        short_rates: np.ndarray,
                        dt: float) -> dict:
        """
        Calibrates the model parameters to observed short rates using Maximum Likelihood Estimation (MLE).
        Args:
            initial_params (VasicekParams): Initial guess for the parameters.
            r0 (float): Initial short rate.
            short_rates (np.ndarray): Observed short rates (not including r0).
            dt (float): Time step between observed rates.
        Returns:
            result (dict): Dictionary containing calibrated parameters and success flag.
        """

        def neg_likelihood(params):
            gamma, barR, sigma = params
            if gamma <= 0 or sigma <= 0:
                return np.inf
            return -self.__likelihood(gamma, barR, sigma, r0, short_rates, dt)

        initial_guess = [initial_params.gamma, initial_params.barR, initial_params.sigma]
        bounds = [(1e-5, None), (None, None), (1e-5, None)]
        result = minimize(neg_likelihood, initial_guess, bounds=bounds)
        calibrated_params = VasicekParams(gamma=result.x[0], barR=result.x[1], sigma=result.x[2])
        return {"result": calibrated_params, 
                "success": result.success}

class VasicekValidator:
    def params_recovery(self, 
                        true_params: VasicekParams, 
                        r0: float,
                        T: float,
                        dt: float, 
                        N: int) -> dict:
        """
        Validates the calibration by checking if the calibrated parameters can recover the true parameters.
        Args:
            true_params (VasicekParams): The true parameters used for simulation.
        """

        # Step1: Simulate short rate paths using true parameters
        simulator = VasicekSimulator().simulate(true_params, r0, T, dt, N, method="transition")
        
        # Step2: Calibrate the model to the simulated short rates
        # storage the calibrated parameters
        estimated_params = []
        initial_params = VasicekParams(gamma=5, barR=0.3, sigma=2)
        for i in range(N):
            rates = simulator[:, i]
            calibrator = VasicekCalibrator()
            res = calibrator.calibrate("mle", 
                                       r0=r0,
                                       initial_params=initial_params,
                                       rates=rates, dt=dt)
            estimated_params.append(res["result"])
        
        # Step3: Compare the calibrated parameters with the true parameters
        dicts = [asdict(p) for p in estimated_params]
        params = np.array([[d['gamma'], d['barR'], d['sigma']] for d in dicts])
        # Calculate bias and std
        bias = np.mean(params, axis=0) - np.array([true_params.gamma, true_params.barR, true_params.sigma])
        std = np.std(params, axis=0)
        # Calculate RMSE, RRMSE
        rmse = np.sqrt(np.mean((params - np.array([true_params.gamma, true_params.barR, true_params.sigma]))**2, axis=0))
        rrmse = rmse / np.array([true_params.gamma, true_params.barR, true_params.sigma])
        return {"bias": bias, "std": std, "rmse": rmse, "rrmse": rrmse}


    def param_recovery_risk_neutral(self, 
                                    true_params: VasicekParams,
                                    r0: float,
                                    maturities: np.ndarray) -> dict:
        """
        Validates the risk-neutral calibration by checking if the calibrated parameters can recover the true parameters.
        Args:
            true_params (VasicekParams): The true parameters used for simulation.
            r0 (float): Initial short rate.
            maturities (np.ndarray): Maturities of the bonds.
        """
        # Step1: Calculate bond prices using true parameters 
        zcb_prices = np.array([VasicekBondPricer.risk_neutral_price(r0, T,
                                                                    true_params.gamma, 
                                                                    true_params.barR, 
                                                                    true_params.sigma) for T in maturities])
        estimated_params = []
        # Step2: Calibrate the model to the bond prices
        for _ in range(100):
            # add random noise to the bond prices
            noise = np.random.normal(0, 0.001, size=zcb_prices.shape)
            zcb_prices += noise
            initial_params = VasicekParams(gamma=5, barR=0.3, sigma=2)
            calibrator = VasicekCalibrator()
            res = calibrator.calibrate("risk_neutral", initial_params=initial_params,
                                        r0=r0,
                                        zcb_prices=zcb_prices, maturities=maturities)
            estimated_params.append(res["result"])
        # Step3: Compare the calibrated parameters with the true parameters
        dicts = [asdict(p) for p in estimated_params]
        params = np.array([[d['gamma'], d['barR'], d['sigma']] for d in dicts])
        # Calculate bias and std
        bias = np.mean(params, axis=0) - np.array([true_params.gamma, 
                                                   true_params.barR, true_params.sigma])
        std = np.std(params, axis=0)
        # Calculate RMSE, RRMSE
        rmse = np.sqrt(np.mean((params - np.array([true_params.gamma, 
                                                   true_params.barR, 
                                                   true_params.sigma]))**2, axis=0))
        rrmse = rmse / np.array([true_params.gamma, true_params.barR, true_params.sigma])
        return {"bias": bias, "std": std, "rmse": rmse, "rrmse": rrmse}



class Vasicek(BaseModel):
    def __init__(self, ):
        self._params = None 
    
    @property
    def params(self) -> VasicekParams:
        """
        Returns the parameters of the Vasicek model.
        Returns:
            VasicekParams: The parameters of the model.
        """
        
        return self._params


    def simulate(self, 
                 r0: float, 
                 T: float, 
                 dt: float, 
                 N: int, 
                 params: VasicekParams | None = None,
                 method: str = "euler") -> np.ndarray: 
        """
        Returns the simulator for the Vasicek model.
        """

        if r0 < 0 or T <= 0 or dt <= 0:
            raise ValueError("Invalid input values. Ensure r0 >= 0, T > 0, and dt > 0.")       
        
        if (params is None) and (self._params is None):
            raise ValueError("Model parameters are not set. Please provide parameters or calibrate the model first.")
        
        if params is None:
            return VasicekSimulator().simulate(self._params, r0, T, dt, N, method)
        else:
            return VasicekSimulator().simulate(params, r0, T, dt, N, method)

    def calibrate(self, 
                  method: str, 
                  *,
                  r0: float | None = None,
                  zcb_prices: np.ndarray | None = None, 
                  maturities: np.ndarray | None = None, 
                  rates: np.ndarray | None = None, 
                  dt: float | None= None,
                  inplace = False) -> dict:
        """
        Calibrates the model parameters to market data.
        Args:
            method (str): Calibration method ("mle" or "risk_neutral").
            initial_params (VasicekParams): Initial guess for the parameters.
            zcb_prices (np.ndarray): Market prices of zero-coupon bonds with face value 1 (required for "risk_neutral").
            maturities (np.ndarray): Corresponding maturities of the bonds (required for "risk_neutral").
            r0 (float): Initial short rate (required for "risk_neutral").
            rates (np.ndarray): Observed short rates (required for "mle").
            dt (float): Time step between observed rates (required for "mle").
            inplace (bool): If True, updates the model parameters in place.
        Returns:
            result (dict): Dictionary containing calibrated parameters and success flag.

        """

        calibrator = VasicekCalibrator()
        result = calibrator.calibrate(method, 
                                           initial_params=VasicekParams(gamma=3, barR=0.3, sigma=0.3),
                                           r0=r0,
                                           zcb_prices=zcb_prices, 
                                           maturities=maturities, 
                                           rates=rates, 
                                           dt=dt)

        if inplace and result["success"]:
            self._params = result["result"]

        return result


