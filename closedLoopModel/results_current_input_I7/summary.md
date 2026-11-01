# Plant driven by the measured currents: fit on the encoder position

Configuration `config_current_input_I7.yaml`; model `config_vcl_fit.yaml` with the first guess of iteration I5 (`results_vcl_study/summary.json`). The measured motor currents (low-passed at 150 Hz, zero phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the start of the window). Windows from 0.3 s before each move to 0.5 s after the end of its command. Identification: MPCLucy_step_1deg, MPCLucy_1degAz_zeroAz_1degPol__zeroPol (18 moves); validation: MPCLucy_2degAz_2degPol_zeroAz_zeroPol, MPCLucy_2degPol_zeroPol_2degAz_zeroAz (8 moves).

Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, `ring`, dampers `c_*` as in vcl_study.py.

## Screening of the parameters

Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the axis) for a change of +/-30 % of each parameter from the first guess; fitted when above 2 %.

| axis | parameter | first guess | cost change | fitted |
|---|---|---|---|---|
| elevation | scale | 0.9268 | 20176.2 % | yes |
| elevation | sigma0 | 620.4 | 134.9 % | yes |
| elevation | Tc | 0.3277 | 8903.3 % | yes |
| elevation | Ts_ratio | 9.988 | 4.8 % | yes |
| elevation | vs | 0.006637 | 16.2 % | yes |
| elevation | d_x_mm | -0.2687 | 37533.5 % | yes |
| elevation | b_{\theta,h} | 0.05 | 2.2 % | yes |

## Fitted parameters

| axis | parameter | first guess (I5) | start (grid) | fitted |
|---|---|---|---|---|
| elevation | scale | 0.9268 | 0.9 | 0.9739 |
| elevation | sigma0 | 620.4 | 620.4 | 1202 |
| elevation | Tc | 0.3277 | 0.3277 | 0.2594 |
| elevation | Ts_ratio | 9.988 | 4 | 6.453 |
| elevation | vs | 0.006637 | 0.006637 | 0.01817 |
| elevation | d_x_mm | -0.2687 | -0.2687 | -0.3254 |
| elevation | b_{\theta,h} | 0.05 | 0.05 | 0.03995 |

## Position fit

NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | first guess (I5) | fitted | adequate |
|---|---|---|---|---|---|
| identification | azimuth | 8 | 96.2 [91, 98] | 96.2 [91, 98] | - |
| identification | elevation | 10 | 94.0 [89, 96] | 95.9 [92, 97] | - |
| validation | azimuth | 4 | 94.5 [91, 98] | 94.5 [91, 98] | yes |
| validation | elevation | 4 | 85.4 [73, 90] | 86.2 [80, 99] | no |
