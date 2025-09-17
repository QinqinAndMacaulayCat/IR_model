
import sys
import os
import numpy as np
from matplotlib import pyplot as plt

from irmodel.model.Vasicek import VasicekValidator, VasicekParams, Vasicek
from irmodel.data.data_fetcher import SampleDataFetcher



def param_recovery_analysis():
    """
    Perform parameter recovery analysis for the Vasicek model.
    1. MLE parameter recovery from simulated short rate paths.
    2. Risk-neutral parameter recovery from simulated bond prices.
    """

    # Part 1: MLE Parameter Recovery
    validator = VasicekValidator()
    true_params = VasicekParams(gamma=3.64, barR=0.3077, sigma=1.813)
    results = validator.params_recovery(true_params, r0=0.03, T=10, dt=1/252, N=10000)

    print("MLE Parameter Recovery Results:")
    print("Bias:", results["bias"])
    print("Standard Deviation:", results["std"])
    print("RMSE:", results["rmse"])
    print("RRMSE:", results["rrmse"])

    # Part 2: Risk-Neutral Parameter Recovery
    print("\nRisk-Neutral Parameter Recovery Results:")
    validator = VasicekValidator()
    true_params = VasicekParams(gamma=3.64, barR=0.3077, sigma=1.813)
    maturities = np.array([0.25, 0.5, 1, 2, 3, 5, 7, 10])
    results = validator.param_recovery_risk_neutral(true_params, r0=0.03, maturities=maturities)
    print("Bias:", results["bias"])
    print("RMSE:", results["rmse"])
    print("RRMSE:", results["rrmse"])



def out_of_sample_mle():
    """
    Calibrate the Vasicek model to historical short rate data using MLE and evaluate out-of-sample performance.
    1. Fetch historical short rate data.
    2. Calibrate Vasicek model using MLE.
    3. Evaluate out-of-sample short rate predictions.
    """

    # Retrieve sample short rate data
    data_fetcher = SampleDataFetcher()
    short_rate_data = data_fetcher.short_rate()

    # Calibrate Vasicek model using MLE
    train_size = int(0.7 * len(short_rate_data))
    vasicek_model = Vasicek()
    result = vasicek_model.calibrate(method='mle', 
                                     r0=short_rate_data['rt'][0],
                                     rates=short_rate_data.loc[:train_size - 1]['rt'].values, 
                                     dt=1/252, 
                                     inplace=True)

    print("Calibrated Vasicek Parameters:")
    print(result)

    # Out-of-sample evaluation
    r0_out = short_rate_data['rt'][train_size - 1]
    rates_out = short_rate_data.iloc[train_size: ]['rt'].values
    T = (short_rate_data.loc[len(short_rate_data) - 1, :]['DATE'] - short_rate_data.loc[train_size - 1, 'DATE']).days / 365
    simulated = vasicek_model.simulate(
            params=result['result'],
            r0=r0_out, 
            T=T, dt=1/512, N=100, 
            method='transition')

    # Match predictions to actual out-of-sample rates
    maturities = ((short_rate_data.loc[train_size:, 'DATE'] - short_rate_data.loc[train_size - 1, 'DATE']).dt.days / 365).values
    indexes = (maturities * 252).astype(int)
    simulated_means = simulated[indexes].mean(axis=1)
    errors = simulated_means - rates_out
    rmse = np.sqrt(np.mean(errors**2))
    rrmse = rmse / np.mean(rates_out)
    print(f"Out-of-sample RMSE: {rmse}, RRMSE: {rrmse}")

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(short_rate_data['DATE'], short_rate_data['rt'], label='Historical Short Rates', color='blue')
    plt.axvline(x=short_rate_data.loc[train_size - 1, 'DATE'], color='red', linestyle='--', label='Train/Test Split')
    plt.plot(short_rate_data.loc[train_size:, 'DATE'], simulated[indexes][:, :10], label='Vasicek Model Prediction', color='green')
    plt.title('Vasicek Model Calibration and Prediction')
    plt.savefig("results/vasicek_calibration.png")


