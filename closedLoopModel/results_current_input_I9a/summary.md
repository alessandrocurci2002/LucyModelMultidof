# Plant driven by the measured currents: fit on the encoder position

Configuration `config_current_input_I9a.yaml`; model `config_vcl_fit.yaml` with the first guess of iteration I5 (`results_vcl_study/summary.json`). The measured motor currents (low-passed at 150 Hz, zero phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the start of the window). Windows from 0.3 s before each move to 0.5 s after the end of its command. Identification: MPCLucy_step_1deg, MPCLucy_1degAz_zeroAz_1degPol__zeroPol (18 moves); validation: MPCLucy_2degAz_2degPol_zeroAz_zeroPol, MPCLucy_2degPol_zeroPol_2degAz_zeroAz (8 moves).

Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, `ring`, dampers `c_*` as in vcl_study.py.

## Screening of the parameters

Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the axis) for a change of +/-30 % of each parameter from the first guess; fitted when above 2 %.

| axis | parameter | first guess | cost change | fitted |
|---|---|---|---|---|
| azimuth | scale | 0.9062 | 7174.0 % | yes |
| azimuth | sigma0 | 202.8 | 29.6 % | yes |
| azimuth | Tc | 0.2279 | 493.2 % | yes |
| azimuth | Ts_ratio | 1.477 | 0.2 % | no |
| azimuth | vs | 0.001 | 0.1 % | no |
| azimuth | K_cable | 1 | 13.0 % | yes |
| azimuth | theta0_cable_deg | 0 | 53.0 % | yes |
| elevation | scale | 0.9739 | 1439.0 % | yes |
| elevation | sigma0 | 1202 | 193.5 % | yes |
| elevation | Tc | 0.2594 | 554.8 % | yes |
| elevation | Ts_ratio | 6.453 | 76.9 % | yes |
| elevation | vs | 0.01817 | 144.9 % | yes |
| elevation | b_{\theta,h} | 0.03995 | 4.5 % | yes |
| elevation | K_cable | 1 | 16.1 % | yes |
| elevation | theta0_cable_deg | 0 | 88.5 % | yes |

## Fitted parameters

| axis | parameter | first guess (I5) | start (grid) | fitted |
|---|---|---|---|---|
| azimuth | scale | 0.9062 | 0.9 | 0.8562 |
| azimuth | sigma0 | 202.8 | 202.8 | 401.5 |
| azimuth | Tc | 0.2279 | 0.2279 | 0.1634 |
| azimuth | K_cable | 1 | 1 | 16.52 |
| azimuth | theta0_cable_deg | 0 | 0.5 | 0.5289 |
| azimuth | Ts_ratio | 1.477 | 1.477 | 2.22 |
| elevation | scale | 0.9739 | 1 | 0.9354 |
| elevation | sigma0 | 1202 | 1202 | 1855 |
| elevation | Tc | 0.2594 | 0.2594 | 0.2422 |
| elevation | Ts_ratio | 6.453 | 6.453 | 3.821 |
| elevation | vs | 0.01817 | 0.01817 | 0.01129 |
| elevation | b_{\theta,h} | 0.03995 | 0.03995 | 0.3041 |
| elevation | K_cable | 1 | 1 | 12.15 |
| elevation | theta0_cable_deg | 0 | 2 | 0.7252 |

## Position fit

NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | first guess (I5) | fitted | adequate |
|---|---|---|---|---|---|
| identification | azimuth | 8 | 95.7 [89, 98] | 95.4 [94, 99] | - |
| identification | elevation | 10 | 74.3 [16, 87] | 96.5 [94, 98] | - |
| validation | azimuth | 4 | 95.7 [95, 96] | 78.7 [77, 80] | no |
| validation | elevation | 4 | 60.0 [51, 84] | 71.5 [61, 85] | no |
