from BaseModel import BaseModel

import numpy as np
from dataclasses import dataclass

@dataclass
class CIRparams:
    gamma: float # Speed of mean reversion
    barR: float # Long-term mean
    alpha: float # Volatility coefficient


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
            self.__simulate_euler(params, r0, T, dt, N)

        elif method == 'transition':
            ...
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






class CIR(BaseModel):

    def __init__(self):
        self._params = None

