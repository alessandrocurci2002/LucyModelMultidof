# Identifiability of the parameters of the current-driven model

Jacobian of the normalised position residuals on the identification moves (1 deg sequences) at the fitted point; fit coordinates: ln p for positive parameters, p for signed ones. Condition number of the Jacobian with columns scaled to unit norm; correlations and standard errors of the estimates from (J^T J)^-1, with the residual variance and the effective number of independent samples (autocorrelation of the residuals). Relative standard error: of p for ln-parameters, absolute for signed ones.

## first run (no cables): azimuth

8 moves, 3936 samples, effective independent samples 27; condition number 14.4; scaled singular values 1.89, 0.513, 0.364, 0.131.

| parameter | fitted value | standard error | correlations |r| > 0.9 |
|---|---|---|---|
| scale | 0.9062 | 47 % | - |
| sigma0 | 202.8 | 587 % | - |
| Tc | 0.2279 | 112 % | - |
| Ts_ratio | 1.477 | 8812 % | - |

## I7 (no cables): elevation

10 moves, 4898 samples, effective independent samples 64; condition number 224; scaled singular values 2.38, 1.02, 0.509, 0.16, 0.0846, 0.0362, 0.0106.

| parameter | fitted value | standard error | correlations |r| > 0.9 |
|---|---|---|---|
| scale | 0.9739 | 24 % | - |
| sigma0 | 1202 | 94 % | - |
| Tc | 0.2594 | 388 % | b_{\theta,h} (-0.99) |
| Ts_ratio | 6.453 | 514 % | - |
| vs | 0.01817 | 380 % | - |
| d_x_mm | -0.3254 | 0.12 (abs.) | - |
| b_{\theta,h} | 0.03995 | 28908 % | Tc (-0.99) |

## I9a (cable spring, one theta_0): azimuth

8 moves, 3936 samples, effective independent samples 35; condition number 17.5; scaled singular values 1.78, 1.15, 0.868, 0.847, 0.179, 0.102.

| parameter | fitted value | standard error | correlations |r| > 0.9 |
|---|---|---|---|
| scale | 0.8562 | 35 % | Tc (+0.90) |
| sigma0 | 401.5 | 2106 % | K_cable (+0.92) |
| Tc | 0.1634 | 173 % | scale (+0.90) |
| Ts_ratio | 2.22 | 47600 % | - |
| K_cable | 16.52 | 737 % | sigma0 (+0.92) |
| theta0_cable_deg | 0.5289 | 0.22 (abs.) | - |

## I9a (cable spring, one theta_0): elevation

10 moves, 4898 samples, effective independent samples 41; condition number 215; scaled singular values 2.44, 1.28, 0.588, 0.229, 0.113, 0.0439, 0.0406, 0.0113.

| parameter | fitted value | standard error | correlations |r| > 0.9 |
|---|---|---|---|
| scale | 0.9354 | 69 % | - |
| sigma0 | 1855 | 237 % | - |
| Tc | 0.2422 | 583 % | b_{\theta,h} (-0.99) |
| Ts_ratio | 3.821 | 1337 % | - |
| vs | 0.01129 | 1273 % | - |
| b_{\theta,h} | 0.3041 | 5384 % | Tc (-0.99) |
| K_cable | 12.15 | 596 % | theta0_cable_deg (-1.00) |
| theta0_cable_deg | 0.7252 | 1.4 (abs.) | K_cable (-1.00) |
