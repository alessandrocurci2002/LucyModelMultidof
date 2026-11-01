# Tracking of the circle with EtherCAT at 1 ms and the ADIS16507 gyro

Configuration `config_dob_ethercat.yaml` on `config_dob_tracking.yaml`. EtherCAT cycle 1 ms, targets linear, feedforward current applied 1 cycle later. Gyro ADIS16507-1 (x, y): 480 Hz second order, 2000 SPS, latency 0.5 ms, ARW 0.13 deg/sqrt(h), bias instability 2.3 deg/h (output noise of the model 0.069 deg/s rms; data sheet 0.068). Metrics on the line of sight from 5 s, one orbit; jitter: rms above each frequency.

## EtherCAT and gyro

| case | E rms [mdeg] | E max [mdeg] | jitter > 1 Hz [urad] | jitter > 10 Hz [urad] | jitter > 100 Hz [urad] | FF current rms az / el [A] | motor current rms az / el [A] |
|---|---|---|---|---|---|---|---|
| current (encoders) | 124.11 | 183.67 | 98.9 | 1.88 | 0.000463 | 0.00 / 0.00 | 2.60 / 1.67 |
| current + gyro observer 20 Hz | 117.66 | 182.98 | 9.66 | 2.63 | 0.142 | 2.86 / 1.35 | 2.85 / 1.35 |
| current + velocity FF | 26.62 | 131.00 | 98.9 | 1.88 | 0.000463 | 0.00 / 0.00 | 2.60 / 1.67 |
| current + velocity FF + gyro observer 5 Hz | 3.90 | 15.20 | 42.5 | 2.92 | 0.0344 | 2.47 / 1.63 | 2.52 / 1.63 |
| current + velocity FF + gyro observer 10 Hz | 2.93 | 10.09 | 18.4 | 2.41 | 0.0695 | 2.62 / 1.55 | 2.62 / 1.55 |
| current + velocity FF + gyro observer 20 Hz | 2.72 | 8.00 | 9.65 | 2.61 | 0.142 | 2.86 / 1.35 | 2.85 / 1.35 |
| current + velocity FF + gyro observer 40 Hz | 2.78 | 8.54 | 10.3 | 3.2 | 0.29 | 3.72 / 1.21 | 3.67 / 1.23 |
| as dob5_ff, gyro without noise | 3.63 | 15.01 | 33.9 | 2.25 | 0.000407 | 2.52 / 1.67 | 2.55 / 1.67 |
| as dob10_ff, gyro without noise | 2.84 | 9.84 | 12 | 1.15 | 0.000164 | 2.55 / 1.66 | 2.55 / 1.66 |
| as dob20_ff, gyro without noise | 2.70 | 7.73 | 5.04 | 0.0112 | 3.14e-06 | 2.55 / 1.66 | 2.54 / 1.66 |
| as dob40_ff, gyro without noise | 2.70 | 7.82 | 4.08 | 0.00822 | 2.46e-06 | 2.55 / 1.66 | 2.54 / 1.66 |
| current + velocity FF, targets held (zoh) | 26.70 | 131.19 | 98.9 | 1.88 | 0.000463 | 0.00 / 0.00 | 2.60 / 1.67 |
| current + velocity FF + gyro observer 20 Hz, targets held (zoh) | 2.90 | 7.99 | 9.65 | 2.61 | 0.142 | 2.86 / 1.35 | 2.85 / 1.35 |

## Continuous loop with an ideal gyro -> EtherCAT and ADIS16507

| case | E rms [mdeg] | E max [mdeg] | jitter > 1 Hz [urad] | jitter > 10 Hz [urad] | jitter > 100 Hz [urad] |
|---|---|---|---|---|---|
| encoders | 124.11 -> 124.11 | 183.67 -> 183.67 | 99.4 -> 98.9 | 1.89 -> 1.88 | 0.00047 -> 0.000463 |
| dob20 | 117.66 -> 117.66 | 182.61 -> 182.98 | 4.42 -> 9.66 | 0.0101 -> 2.63 | 7.13e-06 -> 0.142 |
| encoders_ff | 26.61 -> 26.62 | 130.93 -> 131.00 | 99.4 -> 98.9 | 1.89 -> 1.88 | 0.000471 -> 0.000463 |
| dob20_ff | 2.70 -> 2.72 | 7.79 -> 8.00 | 4.42 -> 9.65 | 0.0101 -> 2.61 | 7.05e-06 -> 0.142 |

## Runs

| case | result | cpu [s] |
|---|---|---|
| dob40_ff_nonoise | completed | 1034 |
| dob20_ff_nonoise | completed | 1041 |
| dob10_ff_nonoise | completed | 1341 |
| dob5_ff_nonoise | completed | 1365 |
| encoders | completed | 243 |
| encoders_ff | completed | 361 |
| dob5_ff | completed | 2934 |
| dob10_ff | completed | 3255 |
| dob20 | completed | 3486 |
| dob20_ff_zoh | completed | 3491 |
| dob20_ff | completed | 3484 |
| encoders_ff_zoh | completed | 923 |
| dob40_ff | completed | 3756 |
