# Plant driven by the measured currents, identified on VClosedLoop/V2

Configuration `config_v2_identification.yaml`, model `config_vcl_fit.yaml` (H6 with the plant of vcl_study.py). The measured motor currents (low-passed at 150 Hz, zero phase) drive the plant without the drives; the goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the window start), windows from 0.3 s before each move to 0.5 s after the end of its command.

Identification: 15 doublet files of V2 (25 azimuth and 24 elevation moves). Excluded: `V2/MPCLucy_doublets_azimut_1deg_15deg_v2` (one move only (-15 -> 0 deg), a repetition of the end of azimut_1deg_15deg). Validation: the 4 step files of VClosedLoop (26 moves, another session).

## 0. Model-free analysis of V2

Rests: mean motor torque K_t i from 3 s after a move to the next move, regressed on angle and direction of the last move (T = K theta + c + h d). Cruise at 10 deg/s: T = F d + K theta + c over all the samples; F is the friction at 10 deg/s, K theta + c the position-dependent torque (torques in nominal K_t units).

| axis | cruise: friction F [N m] | cruise: gradient K [N m/rad] | cruise: offset c [N m] | cruise: rms [N m] | rests: direction term h [N m] | rests: gradient K [N m/rad] | rests: rms [N m] |
|---|---|---|---|---|---|---|---|
| azimuth | 0.374 | 0.160 | +0.049 | 0.510 | 0.201 | 1.121 | 0.097 |
| elevation | 0.384 | 1.289 | +0.074 | 0.109 | 0.305 | 1.624 | 0.155 |

## 1. Rigid surrogate

Each axis turning as a rigid body (inertia r^T M r, same motor, friction and extra torques) against the full H6 model, at the start values, on all the identification windows.

| axis | fit difference, median [points] | max [points] | max position difference / move | wall time rigid [s] | wall time full [s] |
|---|---|---|---|---|---|
| azimuth | 0.02 | 0.59 | 3.52 % | 6 | 17 |
| elevation | 0.03 | 0.26 | 0.32 % | 6 | 14 |

## 2. Structures and cross-validation

Fit on all of V2 (identification median) and leave-one-file-out cross-validation (median of the fits of the held-out moves, each fold refitted from the fit on all of V2); rigid surrogate. Selected: best CV median, or the structure with fewer free parameters within 1 points.

| axis | structure | free parameters | identification median [%] | CV median [%] | CV worst [%] | CV mean (1 - fit)^2 |
|---|---|---|---|---|---|---|
| azimuth | S0 no spring | 5 | 81.5 | 81.1 | -97 | 0.237 |
| azimuth | S1 spring | 7 | 89.0 | 88.8 | -37 | 0.217 |
| azimuth | S1 sigma0 fixed x0.333 | 6 | 88.5 | 88.5 | -46 | 0.221 |
| azimuth | S1 sigma0 fixed x1 | 6 | 81.6 | 81.6 | -5 | 0.146 |
| azimuth | S1 sigma0 fixed x3 **(selected)** | 6 | 89.5 | 89.5 | -23 | 0.120 |
| azimuth | S1 vs fixed x0.333 | 6 | 84.4 | 84.5 | -144 | 0.397 |
| azimuth | S1 vs fixed x1 | 6 | 82.9 | 82.9 | -102 | 0.262 |
| azimuth | S1 vs fixed x3 | 6 | 88.5 | 86.2 | -19 | 0.182 |
| elevation | S0 no spring | 6 | 52.5 | 36.5 | -879 | 4.808 |
| elevation | S1 spring | 8 | 82.9 | 84.1 | -354 | 1.360 |
| elevation | S1 sigma0 fixed x0.333 | 7 | 85.9 | 86.3 | 5 | 0.138 |
| elevation | S1 sigma0 fixed x1 | 7 | 81.8 | 79.6 | -16 | 0.138 |
| elevation | S1 sigma0 fixed x3 | 7 | 82.5 | 82.5 | -25 | 0.207 |
| elevation | S1 vs fixed x0.333 | 7 | 84.2 | 85.8 | -159 | 0.472 |
| elevation | S1 vs fixed x1 | 7 | 87.5 | 86.4 | -161 | 0.470 |
| elevation | S1 vs fixed x3 **(selected)** | 7 | 87.4 | 86.9 | -348 | 1.334 |
| elevation | S1 b_{\theta,h} fixed 0 | 7 | 85.6 | 85.5 | -187 | 0.607 |
| elevation | S1 sigma0, vs, b_{\theta,h} fixed | 5 | 85.1 | 85.4 | -391 | 1.048 |

## 3. Selected structure, full H6 model

Standard errors from the Jacobian of the position residuals with the effective number of independent samples (relative for the positive parameters, absolute for theta0).

| axis | parameter | value | how | standard error | jackknife (CV folds) | rigid surrogate value |
|---|---|---|---|---|---|---|
| azimuth | K_cable | 0.9279 | fitted | 264 % | 20 % | 0.9279 |
| azimuth | Tc | 0.2316 | fitted | 34 % | 12 % | 0.2316 |
| azimuth | Ts_ratio | 7.934 | fitted | 348 % | 11 % | 7.934 |
| azimuth | scale | 0.9644 | fitted | 30 % | 4 % | 0.9644 |
| azimuth | sigma0 | 4513 | fixed (nominal) | - | - | - |
| azimuth | theta0_cable_deg | -2.678 | fitted | 12 deg | 1.4 deg | -2.678 |
| azimuth | vs | 0.005816 | fitted | 116 % | 37 % | 0.005816 |
| elevation | K_cable | 1.661 | fitted | 83 % | 14 % | 1.635 |
| elevation | Tc | 0.4073 | fitted | 10 % | 6 % | 0.3966 |
| elevation | Ts_ratio | 1.139 | fitted | 151 % | 8 % | 1.143 |
| elevation | b_{\theta,h} | 0.0004712 | fitted | 344615 % | 159 % | 0.0004712 |
| elevation | scale | 1.071 | fitted | 38 % | 5 % | 1.053 |
| elevation | sigma0 | 1703 | fitted | 330 % | 97 % | 1674 |
| elevation | theta0_cable_deg | -2.156 | fitted | 3.2 deg | 0.59 deg | -2.136 |
| elevation | vs | 0.015 | fixed (nominal) | - | - | - |

## 4. Position fit, full model

NRMSE fit [%], median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | fit | adequate |
|---|---|---|---|---|
| identification (V2) | azimuth | 25 | 89.6 [-3, 97] | - |
| identification (V2) | elevation | 24 | 87.8 [27, 96] | - |
| validation (VClosedLoop steps) | azimuth | 12 | 85.6 [67, 98] | no |
| validation (VClosedLoop steps) | elevation | 14 | 49.7 [-83, 59] | no |
| validation, theta0 of the session from one other file | azimuth | 36 | 89.5 [76, 97] | no |
| validation, theta0 of the session from one other file | elevation | 42 | 85.4 [77, 94] | no |
| validation, theta0 per sequence | azimuth | 12 | 90.2 [86, 93] | yes |
| validation, theta0 per sequence | elevation | 14 | 86.4 [78, 94] | no |

theta0 of the session: theta0 estimated alone on one validation file (all other parameters from V2) and used on the moves of the other three files; every file in turn (moves counted once per calibration file). theta0 per sequence: estimated on the moves it is then scored on (an unknown initial condition).

| axis | sequence | theta0 [deg] |
|---|---|---|
| azimuth | step_1deg | +2.21 |
| azimuth | 1degAz_zeroAz_1degPol__zeroPol | +0.00 |
| azimuth | 2degAz_2degPol_zeroAz_zeroPol | -0.00 |
| azimuth | 2degPol_zeroPol_2degAz_zeroAz | -0.51 |
| elevation | step_1deg | +2.28 |
| elevation | 1degAz_zeroAz_1degPol__zeroPol | +2.16 |
| elevation | 2degAz_2degPol_zeroAz_zeroPol | +2.86 |
| elevation | 2degPol_zeroPol_2degAz_zeroAz | +2.57 |

theta0 identified on V2: azimuth -2.68 deg, elevation -2.16 deg.

Validation moves:

| move | amplitude [deg] | fit [%] |
|---|---|---|
| step_1deg azimuth 5.76 s | +1.00 | 72.2 |
| step_1deg azimuth 8.04 s | -1.00 | 87.1 |
| step_1deg azimuth 10.37 s | +1.00 | 67.2 |
| step_1deg azimuth 12.01 s | -1.00 | 86.8 |
| step_1deg elevation 15.77 s | +0.57 | 54.9 |
| step_1deg elevation 18.95 s | -1.00 | 30.3 |
| step_1deg elevation 20.16 s | +1.00 | -83.5 |
| step_1deg elevation 21.30 s | -1.00 | 48.7 |
| step_1deg elevation 22.50 s | +1.00 | -64.1 |
| step_1deg elevation 23.32 s | -1.00 | 55.0 |
| 1degAz_zeroAz_1degPol__zeroPol azimuth 3.57 s | +1.00 | 77.5 |
| 1degAz_zeroAz_1degPol__zeroPol azimuth 10.29 s | -1.00 | 96.6 |
| 1degAz_zeroAz_1degPol__zeroPol azimuth 27.94 s | +1.00 | 74.1 |
| 1degAz_zeroAz_1degPol__zeroPol azimuth 33.30 s | -1.00 | 96.6 |
| 1degAz_zeroAz_1degPol__zeroPol elevation 17.00 s | +1.00 | -21.5 |
| 1degAz_zeroAz_1degPol__zeroPol elevation 21.87 s | -1.00 | 54.7 |
| 1degAz_zeroAz_1degPol__zeroPol elevation 37.85 s | +1.00 | -1.0 |
| 1degAz_zeroAz_1degPol__zeroPol elevation 58.75 s | -1.00 | 56.3 |
| 2degAz_2degPol_zeroAz_zeroPol azimuth 18.31 s | +2.00 | 84.5 |
| 2degAz_2degPol_zeroAz_zeroPol azimuth 37.70 s | -2.00 | 97.6 |
| 2degAz_2degPol_zeroAz_zeroPol elevation 23.56 s | +2.00 | 47.1 |
| 2degAz_2degPol_zeroAz_zeroPol elevation 40.79 s | -2.00 | 50.7 |
| 2degPol_zeroPol_2degAz_zeroAz azimuth 19.95 s | +2.00 | 83.9 |
| 2degPol_zeroPol_2degAz_zeroAz azimuth 22.95 s | -2.00 | 97.2 |
| 2degPol_zeroPol_2degAz_zeroAz elevation 8.19 s | +2.00 | 54.3 |
| 2degPol_zeroPol_2degAz_zeroAz elevation 14.17 s | -2.00 | 58.6 |
