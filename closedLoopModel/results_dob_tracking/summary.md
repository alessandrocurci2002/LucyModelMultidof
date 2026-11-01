# Tracking of the circle: current control against the gyro disturbance observer

Configuration `config_dob_tracking.yaml`. Reference of Attorrese (2026): circle of radius 10 m at 5 m, centre 100 m away, 0.2 rad/s. Metrics on the line of sight of the telescope from 5 s (1 orbit). Jitter: rms of the line-of-sight error above the given frequency (1 mdeg = 17.45 urad). Observer: J' = azimuth 2.162 kg m^2, elevation 0.900 kg m^2; K_t' = 0.258 N m/A (true torque scale azimuth 0.96, elevation 1.07).

## Pointing error and jitter

| case | E rms [mdeg] | E mean [mdeg] | E max [mdeg] | miss at the drone, rms [mm] | jitter > 1 Hz [urad] | jitter > 10 Hz [urad] | cross-elevation rms [mdeg] | elevation rms [mdeg] | current az rms [A] | current el rms [A] |
|---|---|---|---|---|---|---|---|---|---|---|
| current: loops on the encoders | 124.11 | 116.69 | 183.67 | 215.9 | 99.4 | 1.89 | 123.82 | 8.48 | 2.60 | 1.67 |
| encoders + gyro observer 5 Hz | 117.66 | 106.71 | 182.61 | 204.7 | 30.5 | 2.03 | 117.51 | 6.03 | 2.55 | 1.67 |
| encoders + gyro observer 10 Hz | 117.66 | 106.71 | 182.61 | 204.7 | 7.93 | 0.586 | 117.50 | 6.04 | 2.55 | 1.66 |
| encoders + gyro observer 20 Hz | 117.66 | 106.71 | 182.61 | 204.7 | 4.42 | 0.0101 | 117.51 | 6.04 | 2.54 | 1.66 |
| encoders + gyro observer 40 Hz | 117.67 | 106.71 | 182.61 | 204.7 | 3.54 | 0.00888 | 117.51 | 6.04 | 2.54 | 1.66 |
| current + velocity feedforward | 26.61 | 14.36 | 130.93 | 46.3 | 99.4 | 1.89 | 26.22 | 4.59 | 2.60 | 1.67 |
| current + velocity feedforward + gyro observer 20 Hz | 2.70 | 1.95 | 7.79 | 4.7 | 4.42 | 0.0101 | 2.70 | 0.05 | 2.54 | 1.66 |
| as above, gyro bandwidth 100 Hz | 2.70 | 1.95 | 7.79 | 4.7 | 4.41 | 0.00957 | 2.70 | 0.05 | 2.54 | 1.66 |
| current + velocity feedforward, alternative azimuth friction | 11.79 | 7.54 | 41.66 | 20.5 | 3.61 | 0.00275 | 10.87 | 4.59 | 1.75 | 1.67 |
| current + velocity feedforward + gyro observer 20 Hz, alternative azimuth friction | 1.68 | 1.55 | 2.97 | 3.0 | 0.407 | 0.00132 | 1.68 | 0.05 | 1.75 | 1.66 |

## Stability of the observer loop

Largest real part of the eigenvalues [1/s] of the linearised loop with the observer on one axis.

| axis | notch on the feedforward | observer bandwidth [Hz] | max real part [1/s] | stable | unstable modes [Hz] |
|---|---|---|---|---|---|
| azimuth | no | 2 | +0.00 | yes | - |
| azimuth | no | 5 | -0.00 | yes | - |
| azimuth | no | 7 | +0.00 | yes | - |
| azimuth | no | 8 | +0.00 | yes | - |
| azimuth | no | 10 | -0.00 | yes | - |
| azimuth | no | 20 | +0.00 | yes | - |
| azimuth | no | 40 | +0.00 | yes | - |
| azimuth | no | 80 | +0.00 | yes | - |
| elevation | yes | 2 | +0.00 | yes | - |
| elevation | yes | 5 | +0.00 | yes | - |
| elevation | yes | 7 | +0.00 | yes | - |
| elevation | yes | 8 | +0.00 | yes | - |
| elevation | yes | 10 | -0.00 | yes | - |
| elevation | yes | 20 | +0.00 | yes | - |
| elevation | yes | 40 | +0.00 | yes | - |
| elevation | yes | 80 | +1.13 | no | 221 |
| elevation | no | 2 | +0.00 | yes | - |
| elevation | no | 5 | +0.00 | yes | - |
| elevation | no | 7 | -0.00 | yes | - |
| elevation | no | 8 | +2.96 | no | 217 |
| elevation | no | 10 | +11.47 | no | 216 |
| elevation | no | 20 | +53.38 | no | 215 |
| elevation | no | 40 | +132.40 | no | 215 |
| elevation | no | 80 | +265.32 | no | 218 |
