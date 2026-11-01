# Resisting torque against the angle and the history (V2)

Motor torque K_t i of the V2 files in the order they were recorded (file modification times). Rests: torque 3 s after the move (start) and at the end of the rest; relaxation T(t) = T_inf + dT exp(-t/tau) from 3 s after the move; integral action: change of torque the velocity integrator of the drive gives from the measured position error, K_t K_pv 2 pi K_iv K_pp int(e dt) (gains of I5).

## Session

| file | start [min] | duration [s] | axis |
|---|---|---|---|
| azimut_1deg_15deg | 0.0 | 100 | azimuth |
| azimut_1deg_15deg_v2 | 28.3 | 100 | azimuth |
| azimut_1deg_minus1deg | 47.2 | 100 | azimuth |
| azimut_3deg_minus3deg | 50.0 | 100 | azimuth |
| azimut_5deg_minus5deg | 52.2 | 100 | azimuth |
| azimut_8deg_minus8deg | 54.3 | 100 | azimuth |
| azimut_01deg_minus01deg | 56.7 | 97 | azimuth |
| azimut_05deg_minus05deg | 59.9 | 99 | azimuth |
| elevation_01deg_minus01deg | 62.7 | 100 | elevation |
| elevation_05deg_minus05deg | 64.8 | 100 | elevation |
| elevation_1deg_minus1deg | 67.0 | 98 | elevation |
| elevation_3deg_minus3deg | 69.3 | 100 | elevation |
| elevation_5deg_minus5deg | 71.7 | 100 | elevation |
| elevation_8deg_minus8deg | 74.0 | 100 | elevation |
| elevation_10deg_minus10deg | 76.2 | 100 | elevation |
| elevation_15deg_minus15deg | 78.6 | 100 | elevation |

## Rests

| file | axis | rest | angle [deg] | last move [deg] | length [s] | T start [N m] | T end [N m] | tau [s] | error start/end [mdeg] | encoder change [mdeg] | integral action dT [N m] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| azimut_1deg_15deg | azimuth | before the first move | +5.00 | - | 2 | +0.171 | +0.170 | - | +0.20 / +0.20 | -0.00 | +0.006 |
| azimut_1deg_15deg | azimuth | after a move | -5.00 | -10.00 | 9 | -0.387 | -0.366 | - | +0.26 / +0.18 | +0.09 | +0.036 |
| azimut_1deg_15deg | azimuth | after a move | +8.00 | +13.00 | 6 | +0.344 | +0.326 | - | +0.24 / +0.24 | -0.00 | +0.028 |
| azimut_1deg_15deg | azimuth | after a move | -8.00 | -16.00 | 6 | -0.322 | -0.306 | - | +0.23 / +0.14 | +0.09 | +0.019 |
| azimut_1deg_15deg | azimuth | after a move | +10.00 | +18.00 | 8 | +0.433 | +0.402 | - | +0.04 / +0.04 | -0.00 | +0.006 |
| azimut_1deg_15deg | azimuth | after a move | -10.00 | -20.00 | 10 | -0.390 | -0.360 | 6 | +0.46 / +0.32 | +0.14 | +0.069 |
| azimut_1deg_15deg | azimuth | after a move | +15.00 | +25.00 | 4 | +0.383 | +0.366 | - | +0.23 / +0.23 | -0.00 | +0.018 |
| azimut_1deg_15deg | azimuth | after a move | -15.00 | -30.00 | 9 | -0.605 | -0.563 | - | +0.39 / +0.16 | +0.23 | +0.038 |
| azimut_1deg_15deg | azimuth | after a move | -0.00 | +15.00 | 4 | +0.366 | +0.350 | - | +0.00 / +0.01 | -0.01 | +0.000 |
| azimut_1deg_15deg_v2 | azimuth | before the first move | -15.00 | - | 36 | -0.612 | -0.586 | - | +0.13 / +0.12 | +0.01 | +0.089 |
| azimut_1deg_15deg_v2 | azimuth | after a move | -0.00 | +15.00 | 59 | +0.426 | +0.347 | 18 | +0.00 / +0.15 | -0.15 | +0.140 |
| azimut_1deg_minus1deg | azimuth | before the first move | +0.00 | - | 9 | -0.100 | -0.100 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_1deg_minus1deg | azimuth | after a move | +1.00 | +1.00 | 35 | +0.202 | +0.192 | 14 | +0.24 / +0.24 | +0.00 | +0.166 |
| azimut_1deg_minus1deg | azimuth | after a move | -1.00 | -2.00 | 32 | -0.216 | -0.190 | 13 | +0.10 / +0.10 | +0.00 | +0.062 |
| azimut_1deg_minus1deg | azimuth | after a move | +0.00 | +1.00 | 14 | +0.194 | +0.187 | 7 | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_3deg_minus3deg | azimuth | before the first move | +0.00 | - | 4 | +0.181 | +0.181 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_3deg_minus3deg | azimuth | after a move | +3.00 | +3.00 | 39 | +0.326 | +0.297 | 16 | +0.05 / +0.05 | +0.00 | +0.035 |
| azimut_3deg_minus3deg | azimuth | after a move | -3.00 | -6.00 | 41 | -0.392 | -0.315 | 10 | +0.30 / +0.30 | +0.00 | +0.241 |
| azimut_3deg_minus3deg | azimuth | after a move | +0.00 | +3.00 | 4 | +0.356 | +0.341 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_5deg_minus5deg | azimuth | before the first move | +0.00 | - | 2 | +0.312 | +0.312 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_5deg_minus5deg | azimuth | after a move | +5.00 | +5.00 | 44 | +0.288 | +0.252 | 16 | +0.19 / +0.19 | +0.00 | +0.163 |
| azimut_5deg_minus5deg | azimuth | after a move | -5.00 | -10.00 | 34 | -0.407 | -0.333 | 9 | +0.15 / +0.15 | +0.00 | +0.103 |
| azimut_5deg_minus5deg | azimuth | after a move | +0.00 | +5.00 | 8 | +0.356 | +0.337 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_8deg_minus8deg | azimuth | after a move | +8.00 | +1.22 | 32 | +0.414 | +0.364 | 9 | +0.24 / +0.24 | +0.00 | +0.146 |
| azimut_8deg_minus8deg | azimuth | after a move | -8.00 | -16.00 | 41 | -0.452 | -0.366 | 9 | +0.11 / +0.11 | +0.00 | +0.086 |
| azimut_8deg_minus8deg | azimuth | after a move | +0.00 | +8.00 | 14 | +0.376 | +0.344 | 6 | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_01deg_minus01deg | azimuth | before the first move | +0.00 | - | 6 | +0.328 | +0.326 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_01deg_minus01deg | azimuth | after a move | +0.10 | +0.10 | 32 | +0.159 | +0.170 | 13 | +0.09 / +0.09 | +0.00 | +0.059 |
| azimut_01deg_minus01deg | azimuth | after a move | -0.10 | -0.20 | 31 | +0.114 | +0.115 | 2 | +0.25 / +0.25 | +0.00 | +0.150 |
| azimut_01deg_minus01deg | azimuth | after a move | +0.00 | +0.10 | 18 | +0.183 | +0.181 | 36 | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_05deg_minus05deg | azimuth | before the first move | +0.00 | - | 6 | +0.180 | +0.180 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| azimut_05deg_minus05deg | azimuth | after a move | +0.50 | +0.50 | 30 | +0.224 | +0.218 | 10 | +0.12 / +0.12 | +0.00 | +0.072 |
| azimut_05deg_minus05deg | azimuth | after a move | -0.50 | -1.00 | 35 | -0.077 | -0.063 | 13 | +0.22 / +0.22 | +0.00 | +0.151 |
| azimut_05deg_minus05deg | azimuth | after a move | +0.00 | +0.50 | 17 | +0.169 | +0.164 | 5 | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_01deg_minus01deg | elevation | after a move | +0.10 | +0.10 | 44 | +0.449 | +0.386 | 17 | +0.20 / +0.09 | +0.11 | +0.097 |
| elevation_01deg_minus01deg | elevation | after a move | -0.10 | -0.20 | 31 | -0.471 | -0.371 | 31 | +0.21 / +0.25 | -0.04 | +0.229 |
| elevation_01deg_minus01deg | elevation | after a move | +0.00 | +0.10 | 14 | +0.393 | +0.353 | 37 | +0.00 / -0.00 | +0.00 | -0.038 |
| elevation_05deg_minus05deg | elevation | before the first move | +0.00 | - | 8 | +0.321 | +0.320 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_05deg_minus05deg | elevation | after a move | +0.50 | +0.50 | 38 | +0.010 | +0.049 | 82 | -2.56 / +0.12 | -2.68 | +0.038 |
| elevation_05deg_minus05deg | elevation | after a move | -0.50 | -1.00 | 35 | +0.043 | +0.119 | 1 | +1.59 / +0.22 | +1.37 | +0.368 |
| elevation_05deg_minus05deg | elevation | after a move | +0.00 | +0.50 | 9 | +0.052 | -0.037 | - | -1.72 / +0.00 | -1.72 | -0.147 |
| elevation_1deg_minus1deg | elevation | before the first move | +0.00 | - | 3 | -0.020 | -0.019 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_1deg_minus1deg | elevation | after a move | +1.00 | +1.00 | 27 | +0.535 | +0.466 | 28 | +0.24 / +0.24 | +0.00 | +0.116 |
| elevation_1deg_minus1deg | elevation | after a move | -1.00 | -2.00 | 41 | -0.389 | -0.243 | 22 | -0.11 / +0.10 | -0.21 | +0.195 |
| elevation_1deg_minus1deg | elevation | after a move | +0.00 | +1.00 | 17 | +0.471 | +0.435 | 42 | +0.01 / -0.00 | +0.02 | -0.007 |
| elevation_3deg_minus3deg | elevation | before the first move | +0.00 | - | 8 | +0.411 | +0.410 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_3deg_minus3deg | elevation | after a move | +3.00 | +3.00 | 33 | +0.609 | +0.493 | 42 | +0.39 / +0.04 | +0.35 | -0.165 |
| elevation_3deg_minus3deg | elevation | after a move | -3.00 | -6.00 | 34 | -0.466 | -0.332 | 28 | +0.29 / +0.30 | -0.00 | +0.344 |
| elevation_3deg_minus3deg | elevation | after a move | +0.00 | +3.00 | 14 | +0.493 | +0.432 | 54 | +0.00 / -0.34 | +0.34 | -0.129 |
| elevation_5deg_minus5deg | elevation | before the first move | +0.00 | - | 5 | +0.359 | +0.357 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_5deg_minus5deg | elevation | after a move | +5.00 | +5.00 | 35 | +0.639 | +0.469 | 31 | +0.08 / +0.19 | -0.10 | -0.086 |
| elevation_5deg_minus5deg | elevation | after a move | -5.00 | -10.00 | 30 | -0.526 | -0.406 | 28 | +0.15 / +0.15 | -0.00 | +0.164 |
| elevation_5deg_minus5deg | elevation | after a move | +0.00 | +5.00 | 19 | +0.517 | +0.441 | 136 | +0.00 / -0.30 | +0.30 | -0.172 |
| elevation_8deg_minus8deg | elevation | before the first move | +0.00 | - | 3 | +0.374 | +0.372 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_8deg_minus8deg | elevation | after a move | +8.00 | +8.00 | 37 | +0.768 | +0.548 | 37 | +0.24 / +0.17 | +0.07 | -0.087 |
| elevation_8deg_minus8deg | elevation | after a move | -8.00 | -16.00 | 26 | -0.554 | -0.419 | 39 | +0.11 / +0.11 | +0.00 | +0.188 |
| elevation_8deg_minus8deg | elevation | after a move | +0.00 | +8.00 | 20 | +0.529 | +0.453 | 163 | +0.00 / -0.22 | +0.22 | -0.170 |
| elevation_10deg_minus10deg | elevation | before the first move | +0.00 | - | 6 | +0.387 | +0.386 | - | +0.00 / +0.00 | +0.00 | +0.000 |
| elevation_10deg_minus10deg | elevation | after a move | +10.00 | +10.00 | 30 | +0.841 | +0.610 | 36 | -0.25 / -0.31 | +0.05 | -0.271 |
| elevation_10deg_minus10deg | elevation | after a move | -10.00 | -20.00 | 36 | -0.624 | -0.431 | 30 | +0.31 / +0.31 | +0.00 | +0.507 |
| elevation_10deg_minus10deg | elevation | after a move | +0.00 | +10.00 | 14 | +0.532 | +0.455 | 500 | +0.00 / -0.34 | +0.34 | -0.134 |
| elevation_15deg_minus15deg | elevation | before the first move | +0.00 | - | 2 | +0.361 | +0.359 | - | +0.00 / +0.00 | +0.00 | -0.000 |
| elevation_15deg_minus15deg | elevation | after a move | +15.00 | +15.00 | 36 | +0.913 | +0.620 | 36 | -0.11 / -0.11 | +0.00 | -0.122 |
| elevation_15deg_minus15deg | elevation | after a move | -15.00 | -30.00 | 34 | -0.763 | -0.496 | 30 | +0.11 / +0.11 | +0.00 | +0.370 |
| elevation_15deg_minus15deg | elevation | after a move | +0.00 | +15.00 | 13 | +0.544 | +0.459 | 406 | +0.00 / -0.34 | +0.34 | -0.119 |

## Cruise, file by file

T = F d + K theta + c over the cruise samples of each file (both directions).

| axis | file | friction F [N m] | gradient K [N m/rad] | offset c [N m] | rms [N m] |
|---|---|---|---|---|---|
| azimuth | azimut_1deg_15deg | 0.282 | 1.256 | +0.036 | 0.127 |
| azimuth | azimut_8deg_minus8deg | 0.359 | 2.797 | +0.082 | 0.133 |
| elevation | elevation_5deg_minus5deg | 0.339 | 1.776 | +0.074 | 0.056 |
| elevation | elevation_8deg_minus8deg | 0.328 | 1.273 | +0.037 | 0.045 |
| elevation | elevation_10deg_minus10deg | 0.343 | 1.255 | +0.045 | 0.043 |
| elevation | elevation_15deg_minus15deg | 0.375 | 1.403 | +0.064 | 0.046 |

Cruise of all files (encoder speed within 3 % of 10 deg/s): azimuth F 0.286 N m, K 1.314 N m/rad, c +0.033 N m; elevation F 0.355 N m, K 1.384 N m/rad, c +0.059 N m.

## The two sessions on the rests

T(3 s) - K theta = c + h d over the rests after moves of 0.9 deg or more, K of the V2 cruise: c is the position-dependent torque at 0, h the friction left at the stop.

| axis | session | rests | c [N m] | h [N m] | rms [N m] |
|---|---|---|---|---|---|
| azimuth | V2 | 22 | +0.018 | 0.237 | 0.098 |
| azimuth | other session | 6 | +0.000 | 0.177 | 0.029 |
| elevation | V2 | 19 | +0.105 | 0.427 | 0.097 |
| elevation | other session | 8 | -0.057 | 0.413 | 0.007 |

## Neutral angle file by file

theta0 fitted alone on the moves of each V2 file, all other parameters as identified on V2 (full H6 model).

| axis | file | theta0 [deg] | fit median [%] | fit median with theta0 of V2 [%] |
|---|---|---|---|---|
| azimuth | azimut_1deg_15deg | -1.82 | 88.7 | 90.1 |
| azimuth | azimut_1deg_15deg_v2 | -5.31 | 88.5 | 81.2 |
| azimuth | azimut_1deg_minus1deg | +0.00 | 94.0 | 82.5 |
| azimuth | azimut_3deg_minus3deg | -1.74 | 94.8 | 94.3 |
| azimuth | azimut_5deg_minus5deg | -1.68 | 91.1 | 89.6 |
| azimuth | azimut_8deg_minus8deg | -0.51 | 88.6 | 86.4 |
| azimuth | azimut_01deg_minus01deg | -7.50 | 37.8 | 26.2 |
| azimuth | azimut_05deg_minus05deg | -5.95 | 77.1 | 90.3 |
| elevation | elevation_01deg_minus01deg | -1.43 | 71.1 | 68.5 |
| elevation | elevation_05deg_minus05deg | -2.21 | 76.5 | 76.0 |
| elevation | elevation_1deg_minus1deg | -3.28 | 89.1 | 84.3 |
| elevation | elevation_3deg_minus3deg | -2.87 | 86.5 | 88.9 |
| elevation | elevation_5deg_minus5deg | -2.56 | 91.7 | 88.3 |
| elevation | elevation_8deg_minus8deg | -2.50 | 92.9 | 90.2 |
| elevation | elevation_10deg_minus10deg | -2.35 | 94.0 | 93.7 |
| elevation | elevation_15deg_minus15deg | -2.53 | 90.6 | 85.4 |

theta0 identified on all of V2: azimuth -2.68 deg, elevation -2.16 deg; other session (VClosedLoop steps): azimuth +2.21, +0.00, -0.00, -0.51; elevation +2.28, +2.16, +2.86, +2.57 deg.
