# Plant driven by the measured currents: fit on the encoder position

Configuration `config_current_input.yaml`; model `config_vcl_fit.yaml` with the first guess of iteration I5 (`results_vcl_study/summary.json`). The measured motor currents (low-passed at 150 Hz, zero phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the start of the window). Windows from 0.3 s before each move to 0.5 s after the end of its command. Identification: MPCLucy_step_1deg, MPCLucy_1degAz_zeroAz_1degPol__zeroPol (18 moves); validation: MPCLucy_2degAz_2degPol_zeroAz_zeroPol, MPCLucy_2degPol_zeroPol_2degAz_zeroAz (8 moves).

Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, `ring`, dampers `c_*` as in vcl_study.py.

## Screening of the parameters

Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the axis) for a change of +/-30 % of each parameter from the first guess; fitted when above 2 %.

| axis | parameter | first guess | cost change | fitted |
|---|---|---|---|---|
| azimuth | scale | 1 | 725.4 % | yes |
| azimuth | sigma0 | 1504 | 7.4 % | yes |
| azimuth | Tc | 0.1726 | 217.0 % | yes |
| azimuth | Ts_ratio | 1.2 | 16.7 % | yes |
| azimuth | vs | 0.001 | 1.0 % | no |
| azimuth | c_{\psi,1} | 0.3696 | 0.0 % | no |
| azimuth | c_{\psi,2} | 6.587 | 0.0 % | no |
| elevation | scale | 1 | 752.2 % | yes |
| elevation | sigma0 | 5826 | 16.5 % | yes |
| elevation | Tc | 0.385 | 846.0 % | yes |
| elevation | Ts_ratio | 1.2 | 94.3 % | yes |
| elevation | vs | 0.001 | 5.1 % | yes |
| elevation | k_b | 1.972e+05 | 1.1 % | no |
| elevation | ring | 1.923 | 0.3 % | no |
| elevation | c_{\theta,1} | 5.04 | 0.0 % | no |
| elevation | c_{\theta,2} | 0.07975 | 0.1 % | no |
| elevation | c_{\theta,3} | 0.8987 | 0.0 % | no |
| elevation | c_{\theta,4} | 24.6 | 0.3 % | no |

## Fitted parameters

| axis | parameter | first guess (I5) | start (grid) | fitted |
|---|---|---|---|---|
| azimuth | scale | 1 | 0.9 | 0.9062 |
| azimuth | sigma0 | 1504 | 1504 | 202.8 |
| azimuth | Tc | 0.1726 | 0.1726 | 0.2279 |
| azimuth | Ts_ratio | 1.2 | 4 | 1.477 |
| elevation | scale | 1 | 1.2 | 0.9374 |
| elevation | sigma0 | 5826 | 5826 | 413.3 |
| elevation | Tc | 0.385 | 0.385 | 0.3959 |
| elevation | Ts_ratio | 1.2 | 4 | 5.196 |
| elevation | vs | 0.001 | 0.001 | 0.006637 |

## Position fit

NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | first guess (I5) | fitted | adequate |
|---|---|---|---|---|---|
| identification | azimuth | 8 | 81.7 [71, 95] | 96.2 [91, 98] | - |
| identification | elevation | 10 | 62.1 [4, 78] | 86.1 [81, 94] | - |
| validation | azimuth | 4 | 79.5 [76, 91] | 94.5 [91, 98] | yes |
| validation | elevation | 4 | 62.9 [52, 78] | 81.0 [58, 95] | no |
