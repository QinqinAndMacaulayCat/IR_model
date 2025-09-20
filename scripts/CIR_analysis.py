"""
This script performs a comprehensive analysis of the CIR model. it includes:

    1. Parameter recovery using synthetic data.
    2. Out-of-sample forecasting and comparison with actual data.
    3. Moving window estimation to observe parameter stability over time.

"""


import numpy as np
from matplotlib import pyplot as plt
from irmodel.model.CIR import CIRValidator, CIRparams, CIRCalibrator
from irmodel.data.data_fetcher import SampleDataFetcher

def parameter_recovery():
    """
    Perform parameter recovery using synthetic data.

    """

    # True parameters
    true_params = CIRparams(gamma=2.2, barR=0.05, alpha=0.1)
    r0 = 0.03
    N = 100
    T = 20
    dt = 1/252
    
    validator = CIRValidator()
    estimated_params = validator.parameter_recovery(true_params=true_params, r0=r0, N=N, T=T, dt=dt)

    print("True parameters:", true_params)
    print("bias", estimated_params['bias'])
    print("std", estimated_params['std'])
    print("rmse", estimated_params['rmse'])
    print("rrmse", estimated_params['rrmse'])
    
def moving_window_estimation():
    """
    Perform moving window estimation to observe parameter stability over time.

    1. Fetch real-world short rate data.
    2. Apply moving window estimation.
    3. Plot parameter estimates over time.

    """
    # Retrieve sample short rate data
    data_fetcher = SampleDataFetcher()
    short_rate_data = data_fetcher.short_rate()

    # Calibrate CIR model using MLE
    train_size = int(0.7 * len(short_rate_data))
    CIR_model = CIRCalibrator()

    params_lst = []
    for window_start in range(0, len(short_rate_data) - train_size, 2):
        result = CIR_model.calibrate(method='mle', 
                                     initial_params=CIRparams(gamma=2, barR=0.03, alpha=0.1),
                                     r0=short_rate_data['rt'][window_start],
                                     short_rates=short_rate_data.loc[:window_start + train_size - 1]['rt'].values, 
                                     dt=1/252)
        if result['success']:
            params_lst.append(result['result'])

    # Convert list of CIRparams to numpy array for easier plotting
    params_array = np.array([[p.gamma, p.barR, p.alpha] for p in params_lst])
    time_points = np.arange(len(params_array))
    param_names = ['gamma', 'barR', 'alpha']
    plt.figure(figsize=(12, 8))
    for i in range(3):
        plt.subplot(3, 1, i + 1)
        plt.scatter(time_points, params_array[:, i], color='blue', label=param_names[i])
        plt.title(f'Moving Window Estimation of {param_names[i]}')
        plt.ylabel(param_names[i])
        plt.grid()
    plt.savefig('./results/CIR_moving_window_estimation.png')


if __name__ == "__main__":

    moving_window_estimation()


