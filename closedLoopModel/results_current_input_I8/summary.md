# Plant driven by the measured currents: fit on the encoder position

Configuration `config_current_input_I8.yaml`; model `config_vcl_fit.yaml` with the first guess of iteration I5 (`results_vcl_study/summary.json`). The measured motor currents (low-passed at 150 Hz, zero phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the start of the window). Windows from 0.3 s before each move to 0.5 s after the end of its command. Identification: MPCLucy_step_1deg, MPCLucy_1degAz_zeroAz_1degPol__zeroPol (18 moves); validation: MPCLucy_2degAz_2degPol_zeroAz_zeroPol, MPCLucy_2degPol_zeroPol_2degAz_zeroAz (8 moves).

Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, `ring`, dampers `c_*` as in vcl_study.py.

## Screening of the parameters

Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the axis) for a change of +/-30 % of each parameter from the first guess; fitted when above 2 %.

| axis | parameter | first guess | cost change | fitted |
|---|---|---|---|---|
| azimuth | scale | 0.9062 | 10501.7 % | yes |
| azimuth | sigma0 | 202.8 | 44.8 % | yes |
| azimuth | Tc | 0.2279 | 695.2 % | yes |
| azimuth | Ts_ratio | 1.477 | 0.3 % | no |
| azimuth | vs | 0.001 | 0.1 % | no |
| azimuth | H_hyst | 0 | 535.6 % | yes |
| azimuth | delta_hyst | 0.001 | 0.0 % | no |
| elevation | scale | 0.9739 | 29.0 % | yes |
| elevation | Tc | 0.2594 | 31.1 % | yes |
| elevation | Ts_ratio | 6.453 | 25.4 % | yes |
| elevation | vs | 0.01817 | 0.0 % | no |
| elevation | b_{\theta,h} | 0.03995 | 0.0 % | no |
| elevation | trig_a | 0 | 0.1 % | no |
| elevation | trig_b | 0 | 0.1 % | no |
| elevation | trig_c | 0 | 0.1 % | no |

## Fitted parameters

| axis | parameter | first guess (I5) | start (grid) | fitted |
|---|---|---|---|---|
| azimuth | scale | 0.9062 | 0.9 | 0.9051 |
| azimuth | sigma0 | 202.8 | 202.8 | 191.4 |
| azimuth | Tc | 0.2279 | 0.2279 | 0.2212 |
| azimuth | H_hyst | 0 | 0 | 0.006057 |
| azimuth | delta_hyst | 0.001 | 0.001 | 0.0003705 |
| azimuth | Ts_ratio | 1.477 | 4 | 2.593 |
| elevation | scale | 0.9739 | 0.9 | 0.9835 |
| elevation | Tc | 0.2594 | 0.2594 | 0.2716 |
| elevation | Ts_ratio | 6.453 | 1.5 | 1.023 |
| elevation | trig_a | 0 | 0 | -18.4 |
| elevation | trig_b | 0 | 0 | -0.01053 |
| elevation | trig_c | 0 | 0 | 0.2227 |

## Position fit

NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | first guess (I5) | fitted | adequate |
|---|---|---|---|---|---|
| identification | azimuth | 8 | 96.2 [91, 98] | 96.1 [91, 98] | - |
| identification | elevation | 10 | -40.0 [-40, -39] | 96.3 [93, 99] | - |
| validation | azimuth | 4 | 94.5 [91, 98] | 94.4 [91, 98] | yes |
| validation | elevation | 4 | -43.1 [-44, -42] | 62.6 [-0, 88] | no |
