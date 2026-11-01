# Plant driven by the measured currents: fit on the encoder position

Configuration `config_current_input_I6.yaml`; model `config_vcl_fit.yaml` with the first guess of iteration I5 (`results_vcl_study/summary.json`). The measured motor currents (low-passed at 150 Hz, zero phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the start of the window). Windows from 0.3 s before each move to 0.5 s after the end of its command. Identification: MPCLucy_step_1deg, MPCLucy_1degAz_zeroAz_1degPol__zeroPol (18 moves); validation: MPCLucy_2degAz_2degPol_zeroAz_zeroPol, MPCLucy_2degPol_zeroPol_2degAz_zeroAz (8 moves).

Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, `ring`, dampers `c_*` as in vcl_study.py.

## Screening of the parameters

Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the axis) for a change of +/-30 % of each parameter from the first guess; fitted when above 2 %.

| axis | parameter | first guess | cost change | fitted |
|---|---|---|---|---|
| elevation | scale | 0.9374 | 4393.2 % | yes |
| elevation | sigma0 | 413.3 | 92.2 % | yes |
| elevation | Tc | 0.3959 | 1781.0 % | yes |
| elevation | Ts_ratio | 5.196 | 0.8 % | no |
| elevation | vs | 0.006637 | 1.3 % | no |
| elevation | d_x_mm | 0 | 9719.7 % | yes |

## Fitted parameters

| axis | parameter | first guess (I5) | start (grid) | fitted |
|---|---|---|---|---|
| elevation | scale | 0.9374 | 0.9 | 0.9268 |
| elevation | sigma0 | 413.3 | 413.3 | 620.4 |
| elevation | Tc | 0.3959 | 0.3959 | 0.3277 |
| elevation | d_x_mm | 0 | 0 | -0.2687 |
| elevation | Ts_ratio | 5.196 | 4 | 9.988 |

## Position fit

NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; adequate when the validation median reaches 90 %.

| set | axis | moves | first guess (I5) | fitted | adequate |
|---|---|---|---|---|---|
| identification | azimuth | 8 | 96.2 [91, 98] | 96.2 [91, 98] | - |
| identification | elevation | 10 | 86.1 [81, 94] | 94.5 [88, 96] | - |
| validation | azimuth | 4 | 94.5 [91, 98] | 94.5 [91, 98] | yes |
| validation | elevation | 4 | 81.0 [58, 95] | 86.8 [75, 91] | no |
