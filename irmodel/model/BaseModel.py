"""
The BaseModel class serves as a foundational class for other models in the application.
It provides common attributes and methods that can be inherited by subclasses:

"""

from abc import ABC, abstractmethod

class BaseModel(ABC):
    
    @abstractmethod
    def simulate(self, r0, dt, T, N, **kwargs):
        """
        Abstract method to simulate a path. Must be implemented by subclasses.
        """
        ...
    @abstractmethod
    def calibrate(self, *kwargs):
        """
        Abstract method to calibrate the model. Must be implemented by subclasses.
        """
        ...
