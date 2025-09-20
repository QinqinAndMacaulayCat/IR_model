import numpy as np
from matplotlib import pyplot as plt

from irmodel.model.Hull_White import HullWhiteCalibrator, HullWhiteFunc, HullWhiteSimulator, f
from irmodel.data.data_fetcher import SampleDataFetcher

 
def find_poly_degree():
    """
    Function to find the optimal polynomial degree for Hull-White calibration.
    """

    degrees = range(1, 10)
    dt = 1/252
    
    fig, ax = plt.subplots(10, 1, figsize=(10, 30))
    
    rates = -np.log(DFs) / Ts
    target_ts = np.arange(0, Ts[-1] + dt, dt)
    for degree in degrees:
        dt = 1/252
        fs, dfs = f().f_and_df(rates, Ts, target_ts, n=degree)
        ax[degree-1].plot(target_ts, fs, label=f'Degree {degree}')
        ax[degree-1].plot(Ts, rates, 'ro', label='Observed Rates')
        ax[degree-1].plot(target_ts, dfs, label=f'Degree {degree} Derivative')

        ax[degree-1].set_title(f'Polynomial Degree {degree}')
        ax[degree-1].legend()

    fig.savefig('./results/hull_white_poly_degree_analysis.png')

   
# Sample usage
fetcher = SampleDataFetcher()
data = fetcher.hull_white_data()
DFs = data['discount_factors'].values
Ts = data['Maturity'].values
caps = data['Cap Price (x100)'].loc[1:].values
cap_strikes = 1 / (1 + data['Swap Rate'].loc[1:].values * 0.25)
notional = 100
        
calibrator = HullWhiteCalibrator()
res = calibrator.calibrate(DFs, caps, cap_strikes, Ts, notional, initial_guess=(0.1, 0.01))
params = res['params']
dt = 1/252
fs, dfs = f().f_and_df(-np.log(DFs) / Ts, Ts, np.arange(0, Ts[-1] + dt, dt), n=6)
gamma, sigma = params.gamma, params.sigma
params.thetas = HullWhiteFunc.thetas(fs, dfs, gamma, sigma, np.arange(0, Ts[-1] + dt, dt))
params.Ts = np.arange(0, Ts[-1] + dt, dt)

print("Calibration success:", res['success'])
print("Calibrated gamma:", res['params'].gamma)
print("Calibrated sigma:", res['params'].sigma)
print("RMSE:", res['rmse'])
print("RRMSE:", res['rrmse'])

simulator = HullWhiteSimulator()
short_rate_paths = simulator.simulate(res['params'], DFs, Ts, N=1000, dt=1/252, seed=42)

# rebuild yield curves
r0 = fs[0]
curves = np.zeros(len(Ts))
for i in range(len(Ts)):
    df = HullWhiteFunc.Z(gamma, sigma, Ts[i], r0, params.thetas, params.Ts)
    curves[i] = -np.log(df) / Ts[i]

# Plot the yield curves
plt.figure(figsize=(10, 6))
plt.plot(Ts, -np.log(DFs) / Ts, 'ro', label='Observed Yield Curve')
plt.plot(Ts, curves, label='Reconstructed Yield Curve')
plt.xlabel('Maturity (Years)')
plt.savefig('./results/hull_white_yield_curve.png')



