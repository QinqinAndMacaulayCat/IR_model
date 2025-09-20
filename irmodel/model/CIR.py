from irmodel.model.BaseModel import BaseModel

import numpy as np
from scipy.stats import ncx2
from scipy.optimize import minimize
from dataclasses import dataclass, asdict

@dataclass
class CIRparams:
    gamma: float # Speed of mean reversion
    barR: float # Long-term mean
    alpha: float # Volatility coefficient

class CIRfunc:
    @staticmethod
    def c(gamma: float,
          alpha: float, 
          dt: float) -> float:
        return 4 * gamma / (alpha * (1 - np.exp(- gamma * dt)))

    @staticmethod
    def nu(gamma: float, 
           barR: float, 
           alpha: float) -> float:
        return 4 * gamma / alpha * barR
    
    @staticmethod
    def lambda_(r_prev: float, 
                gamma: float, 
                alpha: float, 
                dt: float) -> float:
        return CIRfunc.c(gamma, alpha, dt) * r_prev * np.exp(- gamma * dt)


class CIRsimulator:
    
    def simulate(self, 
                 params: CIRparams, 
                 r0: float, 
                 T: float, 
                 dt: float, 
                 N: int,
                 method: str) -> np.ndarray:
        """
        Simulate interest rate paths using the CIR model.
        Parameters:
            params: CIRparams - Model parameters
            r0: float - Initial interest rate
            T: float - Total time horizon
            dt: float - Time step
            N: int - Number of simulation paths
        Returns:
            np.ndarray - Simulated interest rate paths with shape (steps, N) 
        """
        if method == 'euler':
            return self.__simulate_euler(params, r0, T, dt, N)

        elif method == 'transition':
            return self.__simulate_transition(params, r0, T, dt, N)
        else:
            raise ValueError("Wrong method")



    def __simulate_euler(self, 
                         params: CIRparams, 
                         r0: float,
                         T: float,
                         dt: float,
                         N: int) -> np.ndarray:
        """
        Euler-Maruyama method for CIR model simulation.
        """
        steps = int(T / dt) + 1
        rates = np.zeros((steps, N))
        rates[0] = r0
        
        z = np.random.normal(size=(steps - 1, N))

        for i in range(1, steps):
            rates[i] = rates[i - 1] + params.gamma * (params.barR - rates[i - 1]) * dt \
                       + np.sqrt(params.alpha * rates[i - 1]) * np.sqrt(dt) * z[i - 1]
            rates[i] = np.maximum(rates[i], 0)

        return rates


    def __simulate_transition(self,
                              params: CIRparams,
                              r0: float, 
                              T: float,
                              dt: float,
                              N: int) -> np.ndarray:
        """
        Simulate paths using Transiton Density Method  
        """
        steps = int(T / dt) + 1
        rates = np.zeros((steps, N))
        rates[0] = r0
        nu = CIRfunc.nu(params.gamma, params.barR, params.alpha) 
        c = CIRfunc.c(params.gamma, params.alpha, dt)

        for i in range(1, steps):
            lambda_ = c * rates[i - 1] * np.exp(-params.gamma * dt)
            z = np.random.noncentral_chisquare(nu, lambda_, size=N)
            rates[i] = z / c
            rates[i] = np.maximum(rates[i], 0)  # Ensure non-negativity

        return rates


class CIRCalibrator:

    def calibrate(self, method: str,
                  initial_params: CIRparams,
                  r0: float,
                  dt: float,
                  short_rates: np.ndarray
                  ) -> dict:
        """
        Calibrate CIR model parameters to market data.
        """
        if method == 'mle':
            return self.__calibrate_mle(initial_params, r0, dt, short_rates)
        elif method == 'risk_neutral':
            ...
        else:
            raise ValueError("Wrong calibration method")


    def __likelihood(self, 
                     gamma: float,
                     barR: float,
                     alpha: float,
                     r0: float,
                     dt: float,
                     short_rates: np.ndarray) -> float:

                      
        """
        Compute the likelihood of observed data given model parameters.
        """       

        if gamma <= 0 or barR <= 0 or alpha<= 0:
            return -np.inf

        feller_gap = alpha - 2.0 * gamma * barR
        penalty = 0.0
        if feller_gap > 0:
            penalty += 1e6 * feller_gap * feller_gap 

        eps = 1e-10 # in case of zero
        c = CIRfunc.c(gamma, alpha, dt)       
        nu = CIRfunc.nu(gamma, barR, alpha)
        lambdas = np.zeros(len(short_rates))
        lambdas[0] = CIRfunc.lambda_(r0, gamma, alpha, dt)
        for i in range(1, len(short_rates)):
            r_prev = short_rates[i - 1]
            lambdas[i] = CIRfunc.lambda_(r_prev, gamma, alpha, dt)

        x = np.maximum(c * short_rates, eps)
        lam = np.maximum(lambdas, eps)

        logpdfs = np.log(c) + ncx2.logpdf(x, df=nu, nc=lam)
        ll = np.sum(logpdfs)

        return -ll + penalty

    def __calibrate_mle(self, 
                        initial_params: CIRparams,
                        r0: float,
                        dt: float,
                        short_rates: np.ndarray) -> dict:
        """
        Calibrate parameters using Maximum Likelihood Estimation (MLE).

        """

        def likelihood(params):
            gamma, barR, alpha = params
            return self.__likelihood(gamma, barR, alpha,
                                     r0, dt, short_rates)
        
        initial_guess = [initial_params.gamma, initial_params.barR, initial_params.alpha]
        bounds = [(1e-5, None), (1e-5, None), (1e-5, None)]
        result = minimize(likelihood,
                          initial_guess,
                          bounds=bounds)

        params = CIRparams(gamma=result.x[0], 
                           barR=result.x[1], 
                           alpha=result.x[2])
        return {"result": params, "success": result.success}


class CIRValidator:

    def parameter_recovery(self,
                           true_params: CIRparams,
                           r0: float,
                           T: float,
                           dt: float,
                           N: int) -> dict:

        """
        Perform parameter recovery using synthetic data.
        1. Generate synthetic data using known parameters.
        2. Estimate parameters from the synthetic data.
        3. Compare estimated parameters with true parameters.
        """
        # Step1: Simulate short rate paths using true parameters
        simulator = CIRsimulator().simulate(true_params, r0, T, dt, N, method="transition")
        
        # Step2: Calibrate the model to the simulated short rates
        # storage the calibrated parameters
        estimated_params = []
        initial_params = CIRparams(gamma=5, barR=0.3, alpha=0.2)
        calibrator = CIRCalibrator()
        for i in range(N):
            rates = simulator[:, i]
            res = calibrator.calibrate(method="mle", 
                                       r0=r0,
                                       initial_params=initial_params,
                                       short_rates=rates, dt=dt)
            if res["success"]:
                estimated_params.append(res["result"])
        
        # Step3: Compare the calibrated parameters with the true parameters
        dicts = [asdict(p) for p in estimated_params]
        params = np.array([[d['gamma'], d['barR'], d['alpha']] for d in dicts])
        # Calculate bias and std
        bias = np.mean(params, axis=0) - np.array([true_params.gamma, true_params.barR, true_params.alpha])
        std = np.std(params, axis=0)
        # Calculate RMSE, RRMSE
        rmse = np.sqrt(np.mean((params - np.array([true_params.gamma, true_params.barR, true_params.alpha]))**2, axis=0))
        rrmse = rmse / np.array([true_params.gamma, true_params.barR, true_params.alpha])
        return {"bias": bias, "std": std, "rmse": rmse, "rrmse": rrmse}




class CIR(BaseModel):

    def __init__(self):
        self._params = None

    @property
    def params(self) -> CIRparams | None:
        return self._params

    def simulate(self,
                 r0: float,
                 T: float,
                 dt: float,
                 N: int,
                 params: CIRparams | None = None,
                 method: str = 'transition') -> np.ndarray:
        """
        Simulate interest rate paths using the CIR model.
        Parameters:
            r0: float - Initial interest rate
            T: float - Total time horizon
            dt: float - Time step
            N: int - Number of simulation paths
            method: str - Simulation method ('euler' or 'transition')
        Returns:
            np.ndarray - Simulated interest rate paths with shape (steps, N) 
        """
        if (self._params is None) and (params is None):
            raise ValueError("Model parameters not set. Please calibrate the model first.")

        if params is not None:
            simulator = CIRsimulator()
            return simulator.simulate(params, r0, T, dt, N, method)
        
        if self._params is not None:
            simulator = CIRsimulator()
            return simulator.simulate(self._params, r0, T, dt, N, method)

    def calibrate(self, 
                  method: str,
                  initial_params: CIRparams,
                  r0: float,
                  dt: float,
                  short_rates: np.ndarray,
                  inplace: bool = True
                  ) -> dict:
        """
        Calibrate CIR model parameters to market data.
        Parameters:
            method: str - Calibration method ('mle' or 'risk_neutral')
            initial_params: CIRparams - Initial guess for model parameters
            r0: float - Initial interest rate
            dt: float - Time step
            short_rates: np.ndarray - Observed short rate data
            inplace: bool - If True, update model parameters with calibrated values
        Returns:
            dict - Calibration results including optimized parameters and success flag
        """
        calibrator = CIRCalibrator()
        result = calibrator.calibrate(method, initial_params, r0, dt, short_rates)
        if result["success"] and inplace:
            self._params = result["result"]
        return result



