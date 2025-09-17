
# Vasicek

## Sensitivity Analysis

### Parameter Recovery

We examine the robustness of the calibration results.

1. Simulate short rates and calibrate the model using MLE. 

The result shows that the sigma parameter is relatively stable, while the gamma and bar_r are sensitive to short rate changes.

2. Add noise to ZCB prices and calibrate the model using risk-neutral calibration.

The result shows that the sigma parameter is sensitive to noise, which indicates we should be cautious when using ZCB prices for calibration. We may consider calibrating volatility separately rather than using ZCB prices and risk-neutral calibration.

```plaintext
MLE Parameter Recovery Results:
Bias: [ 0.48054904 -0.01506855  0.00059596]
Standard Deviation: [0.94724384 0.14930257 0.02553081]
RMSE: [1.06216678 0.15006105 0.02553776]
RRMSE: [0.29180406 0.48768622 0.01408591]

Risk-Neutral Parameter Recovery Results:
Bias: [ 0.19455133 -0.06846515 -0.6786019 ]
RMSE: [0.59719612 0.07656089 0.9596563 ]
RRMSE: [0.16406487 0.24881669 0.52931952]
```

### Out-of-Sample Simulation 

We simulate short rates using the calibrated parameters and compare the simulated short rates with the actual short rates.

From the result, we can see the forecast performance is not very good, which indicates the Vasicek model may not be suitable for modeling short rates. This is consistent with the findings in the literature.


![Out-of-Sample Simulation](../scripts/results/vasicek_calibration.png)
