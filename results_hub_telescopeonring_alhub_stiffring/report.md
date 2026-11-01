# Multi-DOF model of the alt-azimuth gimbal: report

> **Realism note (section 0 of the specification).** Shaft geometries,
> materials, joint stiffnesses, rotor and fork inertias and damping values
> are currently **placeholders**. Every frequency, transfer function and figure of
> this report is a **result from placeholder parameters** and must not be presented as a
> characterisation of the real system.


## Reference frames (De Donatis convention)

- N inertial: origin O_N (origin of the Onshape assembly, base of the azimuth shaft), unit vectors x_hat, y_hat, z_hat;
- body 1 (azimuth): origin o_1 on the azimuth axis, dist(O_N, o_1) = 59.695 mm along z_hat; unit vectors i_hat, j_hat, k_hat = z_hat;
- body 2 (polar): origin o_2 = centre of rotation of the polar shaft, dist(O_N, o_2) = 548.865 mm along z_hat; unit vectors K_hat (pointing, towards the optics), J_hat = j_hat (polar axis), I_hat = J_hat x K_hat;
- the origins of the bodies lie on the rotation axes, not at the centres of mass.

## Parameters used

### Materials

| material | rho [kg/m^3] | G [GPa] | placeholder |
|---|---|---|---|
| steel | 8000 | 77 | **yes** |
| aluminium | 2810 | 26.9 | **yes** |

### Azimuth axis

| element | mode | L [mm] | J [kg m^2] | k [N m/rad] | m [kg] | I_t [kg m^2] | placeholder |
|---|---|---|---|---|---|---|---|
| 1 | geometry | 29 | 1.493e-06 | 1.708e+04 | 0.04665 | 4.015e-06 | no |
| 2 | geometry | 20 | 0.0003457 | 8.278e+06 | 0.4463 | 0.0001877 | no |
| 3 | geometry | 22.59 | 1.163e-06 | 2.193e+04 | 0.03634 | 2.127e-06 | no |
| 4 | geometry | 10 | 2.149e-05 | 2.068e+06 | 0.1039 | 1.161e-05 | no |
| 5 | geometry | 10 | 2.149e-05 | 2.068e+06 | 0.1039 | 1.161e-05 | no |
| 6 | geometry | 30.88 | 1.589e-06 | 1.604e+04 | 0.04967 | 4.742e-06 | no |
| 7 | geometry | 35.53 | 1.09e-05 | 8.311e+04 | 0.1395 | 2.013e-05 | no |
| 8 | geometry | 17.5 | 2.063e-05 | 6.482e+05 | 0.1347 | 1.375e-05 | no |

| node | inertia | J [kg m^2] | source | placeholder |
|---|---|---|---|---|
| 9 | load (fork + polar shafts) | 1.252 | Onshape, Tavola Rotante [DEF] <1>, without mate connector | no |

Breakdown by layer: `results/tables/element_breakdown_azimut.csv`.

### Polar axis

| element | mode | L [mm] | J [kg m^2] | k [N m/rad] | m [kg] | I_t [kg m^2] | placeholder |
|---|---|---|---|---|---|---|---|
| 1 | geometry | 58.29 | 1.788e-05 | 5.066e+04 | 0.2289 | 7.375e-05 | no |
| 2 | hub | 214 | 0 | 1.506e+05 | 0 | 0 | no |
| 3 | telescope | 1 | 0 | -- | 0 | 0 | no |
| 4 | hub | 214 | 0 | 1.506e+05 | 0 | 0 | no |
| 5 | geometry | 47.9 | 2.466e-06 | 1.034e+04 | 0.07705 | 1.596e-05 | no |
| 6 | geometry | 10.5 | 2.263e-05 | 1.976e+06 | 0.1093 | 1.232e-05 | no |
| 7 | geometry | 10.5 | 2.263e-05 | 1.976e+06 | 0.1093 | 1.232e-05 | no |
| 8 | geometry | 38.74 | 1.994e-06 | 1.279e+04 | 0.06231 | 8.79e-06 | no |
| 9 | geometry | 27 | 0.001263 | 1.668e+07 | 1.309 | 0.0007111 | no |
| 10 | geometry | 12.08 | 6.22e-07 | 4.1e+04 | 0.01944 | 5.475e-07 | no |

| node | inertia | J [kg m^2] | source | placeholder |
|---|---|---|---|---|

Breakdown by layer: `results/tables/element_breakdown_polare.csv`.

### Telescope (payload, Q3)

- source: config.yaml, frame `cad_cg` (Onshape, Telescopio LAMERA <1>, without mate connector); subtraction of the polar parts: no
- resulting values, about o_2 in the frame of body 2 (I_hat, J_hat, K_hat):
  - mass m_p = 15.3079 kg, eccentricity d_2 = (d_I, d_J, d_K) = (-5.20, -0.00, +73.12) mm (K_hat = pointing)
  - I_II = 0.857821, I_JJ = 0.866208, I_KK = 0.256201 kg m^2; I_IJ = 7.126e-07, I_IK = -7.508e-03, I_JK = -1.139e-06 kg m^2
- subtracted: mass 0.0000 kg, sum J = 0.000e+00 kg m^2
- J^ax = 5.3312e-02 kg m^2 (J^ax = 8.3476e-04 kg m^2 from the polar elements is already in the azimuth load: not added)

## Reduction

Stiffnesses of the reduced chain: k_psi = [ 9547.55077 13090.75942] N m/rad, k_theta = [37907.56488  9630.88359  9694.18968  9800.     ] N m/rad.
Damping: beta_psi = 4.2318e-06 s, beta_theta = 3.6376e-05 s.

### Validation (isolated axes, theta0 = 0)

| axis | mass | reduced mass | mode | f full [Hz] | f reduced [Hz] | error | 3 dB limit f [Hz] |
|---|---|---|---|---|---|---|---|
| azimuth | consistent | guyan | 1 | 0 | 0 | +0.000 % | 659.5 |
| azimuth | consistent | guyan | 2 | 752.185 | 834.269 | +10.913 % |  |
| azimuth | consistent | guyan | 3 | 4595.13 | 6167.57 | +34.220 % |  |
| polar | consistent | guyan | 1 | 0 | 0.148206 | +0.000 % | 319.2 |
| polar | consistent | guyan | 2 | 87.5063 | 87.5192 | +0.015 % |  |
| polar | consistent | guyan | 3 | 336.31 | 354.342 | +5.362 % |  |
| polar | consistent | guyan | 4 | 3506.35 | 3810.05 | +8.662 % |  |
| polar | consistent | guyan | 5 | 10718.9 | 11070.2 | +3.277 % |  |
| azimuth | consistent | guyan_lumped | 1 | 0 | 0 | +0.000 % | 679.4 |
| azimuth | consistent | guyan_lumped | 2 | 752.185 | 811.21 | +7.847 % |  |
| azimuth | consistent | guyan_lumped | 3 | 4595.13 | 1818.41 | -60.427 % |  |
| polar | consistent | guyan_lumped | 1 | 0 | 0.148206 | +0.000 % | 348.1 |
| polar | consistent | guyan_lumped | 2 | 87.5063 | 87.5189 | +0.014 % |  |
| polar | consistent | guyan_lumped | 3 | 336.31 | 346.641 | +3.072 % |  |
| polar | consistent | guyan_lumped | 4 | 3506.35 | 1246.46 | -64.451 % |  |
| polar | consistent | guyan_lumped | 5 | 10718.9 | 9264.02 | -13.573 % |  |
| azimuth | lumped | guyan | 1 | 0 | 0 | +0.000 % | 659.5 |
| azimuth | lumped | guyan | 2 | 752.066 | 834.085 | +10.906 % |  |
| azimuth | lumped | guyan | 3 | 4572.95 | 6108.64 | +33.582 % |  |
| polar | lumped | guyan | 1 | 0.148206 | 0.148206 | +0.000 % | 319.2 |
| polar | lumped | guyan | 2 | 87.5063 | 87.5192 | +0.015 % |  |
| polar | lumped | guyan | 3 | 336.286 | 354.315 | +5.361 % |  |
| polar | lumped | guyan | 4 | 3480.57 | 3781.71 | +8.652 % |  |
| polar | lumped | guyan | 5 | 9988.59 | 10051.4 | +0.628 % |  |
| azimuth | lumped | guyan_lumped | 1 | 0 | 0 | +0.000 % | 679.4 |
| azimuth | lumped | guyan_lumped | 2 | 752.066 | 811.21 | +7.864 % |  |
| azimuth | lumped | guyan_lumped | 3 | 4572.95 | 1818.41 | -60.235 % |  |
| polar | lumped | guyan_lumped | 1 | 0.148206 | 0.148206 | +0.000 % | 348.1 |
| polar | lumped | guyan_lumped | 2 | 87.5063 | 87.5189 | +0.014 % |  |
| polar | lumped | guyan_lumped | 3 | 336.286 | 346.641 | +3.079 % |  |
| polar | lumped | guyan_lumped | 4 | 3480.57 | 1246.46 | -64.188 % |  |
| polar | lumped | guyan_lumped | 5 | 9988.59 | 9264.02 | -7.254 % |  |

## Operating points, frequencies and zeros

| mass | form | theta0 [deg] | Omega0 [deg/s] | tau_theta0 [N m] | first frequencies [Hz] | unstable pole |
|---|---|---|---|---|---|---|
| consistent | guyan | -45 | 0 | 7.2105 | 0.48378, 87.519, 354.34, 834.28, 3810.1 | no |
| consistent | guyan | -45 | 10 | 7.2013 | 0.48397, 87.519, 354.34, 834.28, 3810.1 | no |
| consistent | guyan | -15 | 0 | 10.401 | 0.31812, 87.519, 354.34, 834.27, 3810.1 | no |
| consistent | guyan | -15 | 10 | 10.397 | 0.31888, 87.519, 354.34, 834.27, 3810.1 | no |
| consistent | guyan | 0 | 0 | 10.977 | 0.14821, 87.519, 354.34, 834.27, 3810.1 | no |
| consistent | guyan | 0 | 10 | 10.978 | 0.14994, 87.519, 354.34, 834.27, 3810.1 | no |
| consistent | guyan | 15 | 0 | 10.805 | 87.519, 354.34, 834.27, 3810.1, 6167.6 | 1.523 rad/s |
| consistent | guyan | 15 | 10 | 10.81 | 87.519, 354.34, 834.27, 3810.1, 6167.6 | 1.517 rad/s |
| consistent | guyan | 45 | 0 | 8.3139 | 87.519, 354.34, 834.28, 3810.1, 6167.6 | 2.831 rad/s |
| consistent | guyan | 45 | 10 | 8.3231 | 87.519, 354.34, 834.28, 3810.1, 6167.6 | 2.83 rad/s |
| consistent | guyan_lumped | -45 | 0 | 7.2105 | 0.48378, 87.519, 346.64, 811.22, 1246.5 | no |
| consistent | guyan_lumped | -45 | 10 | 7.2013 | 0.48397, 87.519, 346.64, 811.22, 1246.5 | no |
| consistent | guyan_lumped | -15 | 0 | 10.401 | 0.31812, 87.519, 346.64, 811.21, 1246.5 | no |
| consistent | guyan_lumped | -15 | 10 | 10.397 | 0.31888, 87.519, 346.64, 811.21, 1246.5 | no |
| consistent | guyan_lumped | 0 | 0 | 10.977 | 0.14821, 87.519, 346.64, 811.21, 1246.5 | no |
| consistent | guyan_lumped | 0 | 10 | 10.978 | 0.14994, 87.519, 346.64, 811.21, 1246.5 | no |
| consistent | guyan_lumped | 15 | 0 | 10.805 | 87.519, 346.64, 811.21, 1246.5, 1818.4 | 1.523 rad/s |
| consistent | guyan_lumped | 15 | 10 | 10.81 | 87.519, 346.64, 811.21, 1246.5, 1818.4 | 1.517 rad/s |
| consistent | guyan_lumped | 45 | 0 | 8.3139 | 87.519, 346.64, 811.22, 1246.5, 1818.4 | 2.831 rad/s |
| consistent | guyan_lumped | 45 | 10 | 8.3231 | 87.519, 346.64, 811.22, 1246.5, 1818.4 | 2.83 rad/s |
| lumped | guyan | -45 | 0 | 7.2105 | 0.48378, 87.519, 354.32, 834.1, 3781.7 | no |
| lumped | guyan | -45 | 10 | 7.2013 | 0.48397, 87.519, 354.32, 834.1, 3781.7 | no |
| lumped | guyan | -15 | 0 | 10.401 | 0.31812, 87.519, 354.32, 834.09, 3781.7 | no |
| lumped | guyan | -15 | 10 | 10.397 | 0.31888, 87.519, 354.32, 834.09, 3781.7 | no |
| lumped | guyan | 0 | 0 | 10.977 | 0.14821, 87.519, 354.32, 834.08, 3781.7 | no |
| lumped | guyan | 0 | 10 | 10.978 | 0.14994, 87.519, 354.32, 834.08, 3781.7 | no |
| lumped | guyan | 15 | 0 | 10.805 | 87.519, 354.32, 834.09, 3781.7, 6108.6 | 1.523 rad/s |
| lumped | guyan | 15 | 10 | 10.81 | 87.519, 354.32, 834.09, 3781.7, 6108.6 | 1.517 rad/s |
| lumped | guyan | 45 | 0 | 8.3139 | 87.519, 354.32, 834.1, 3781.7, 6108.6 | 2.831 rad/s |
| lumped | guyan | 45 | 10 | 8.3231 | 87.519, 354.32, 834.1, 3781.7, 6108.6 | 2.83 rad/s |
| lumped | guyan_lumped | -45 | 0 | 7.2105 | 0.48378, 87.519, 346.64, 811.22, 1246.5 | no |
| lumped | guyan_lumped | -45 | 10 | 7.2013 | 0.48397, 87.519, 346.64, 811.22, 1246.5 | no |
| lumped | guyan_lumped | -15 | 0 | 10.401 | 0.31812, 87.519, 346.64, 811.21, 1246.5 | no |
| lumped | guyan_lumped | -15 | 10 | 10.397 | 0.31888, 87.519, 346.64, 811.21, 1246.5 | no |
| lumped | guyan_lumped | 0 | 0 | 10.977 | 0.14821, 87.519, 346.64, 811.21, 1246.5 | no |
| lumped | guyan_lumped | 0 | 10 | 10.978 | 0.14994, 87.519, 346.64, 811.21, 1246.5 | no |
| lumped | guyan_lumped | 15 | 0 | 10.805 | 87.519, 346.64, 811.21, 1246.5, 1818.4 | 1.523 rad/s |
| lumped | guyan_lumped | 15 | 10 | 10.81 | 87.519, 346.64, 811.21, 1246.5, 1818.4 | 1.517 rad/s |
| lumped | guyan_lumped | 45 | 0 | 8.3139 | 87.519, 346.64, 811.22, 1246.5, 1818.4 | 2.831 rad/s |
| lumped | guyan_lumped | 45 | 10 | 8.3231 | 87.519, 346.64, 811.22, 1246.5, 1818.4 | 2.83 rad/s |

Invariant zeros of the torque -> encoder matrix (consistent/guyan, theta0 = -45 deg, Omega0 = 0), frequencies in Hz:

| TF | zeros |
|---|---|
| psi_e/tau_psi | 0.4838, 13.316, 87.519, 354.34, 1642.7, 1717.7 |
| psi_e/tau_theta | 0.0005, 577.35, 1642.8, 1717.8, 4218.6, 4375.4 |
| theta_e/tau_psi | 0, 346.94, 1886, 3809.7, 4116.4, 4375.3 |
| theta_e/tau_theta | 0, 16.936, 577.35, 834.28, 4116.4, 4218.6 |

Full tables: `results/tables/poles_zeros.csv`; state matrices in `results/tf/`.

## Gawronski modal truncation (default reduction)

> Every number in this section is a **result from placeholder parameters**.

Modal basis from the FULL FEM of each axis with the payload frozen at theta_ref = 0 deg (k_g = 0); rigid mode 1 not normalised (eta_0 = angle, modal mass J_tot), flexible modes mass-normalised. Band 2000 Hz, cut-off 4000 Hz, hankel_min_ratio = None, n_flex = {'azimuth': None, 'polar': None}, static correction: none.

### Azimuth axis (J_tot = 2.1637 kg m^2, beta = 4.2318e-06 s)

| mode | f [Hz] | zeta | c_e | c_m | c_L | Hankel [rad/(N m)] | kept | reason |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 1 | 1 | 1 | inf | **yes** | rigid |
| 1 | 752.185 | 0.01 | 52.33 | 33.54 | -0.009113 | 1.964e-03 | **yes** | cut-off |
| 2 | 4595.13 | 0.0611 | -12.59 | 147 | -0.001177 | 9.087e-06 | **yes** | validation |
| 3 | 21343.7 | 0.284 | 0.1758 | -12.45 | -0.00167 | 1.072e-10 | no |  |
| 4 | 29518.9 | 0.392 | 1416 | 0.09145 | 1.104e-07 | 2.398e-09 | no |  |
| 5 | 47898.7 | 0.637 | 0.0255 | 3.527 | -0.001505 | 3.897e-13 | no |  |
| 6 | 83337 | 1.11 | 10.06 | 1.326 | 1.839e-05 | 1.098e-11 | no |  |
| 7 | 85155.5 | 1.13 | 64.66 | -1.388 | -2.712e-06 | 6.922e-11 | no |  |

### Polar axis (J_tot = 0.89979 kg m^2, beta = 3.6376e-05 s)

| mode | f [Hz] | zeta | c_e | c_m | c_P | Hankel [rad/(N m)] | kept | reason |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 1 | 1 | 1 | inf | **yes** | rigid |
| 1 | 87.5063 | 0.01 | 5.338 | 5.574 | -0.2075 | 2.461e-03 | **yes** | cut-off |
| 2 | 336.31 | 0.0384 | -1.175 | 15.18 | 0.002981 | 2.599e-05 | **yes** | cut-off |
| 3 | 3506.35 | 0.401 | -0.1066 | 145.5 | 2.14e-06 | 1.993e-08 | **yes** | cut-off |
| 4 | 10718.9 | 1.22 | 333.9 | 0.007906 | 3.16e-07 | 1.188e-10 | no |  |
| 5 | 39994.3 | 4.57 | -334.9 | 0.0006959 | 6.051e-09 | 2.018e-13 | no |  |
| 6 | 63122.7 | 7.21 | -4.055e-05 | 0.542 | 1.429e-11 | 4.843e-18 | no |  |
| 7 | 68636.4 | 7.84 | -0.00374 | -16.91 | 1.483e-09 | 1.084e-14 | no |  |

### Spectrum of the two axes (rigid modes first)

| # | f [Hz] | axis | mode | kept |
|---|---|---|---|---|
| 1 | 0 | azimuth | rigid | yes |
| 2 | 0 | polar | rigid | yes |
| 3 | 87.5063 | polar | flexible 1 | yes |
| 4 | 336.31 | polar | flexible 2 | yes |
| 5 | 752.185 | azimuth | flexible 1 | yes |
| 6 | 3506.35 | polar | flexible 3 | yes |
| 7 | 4595.13 | azimuth | flexible 2 | yes |
| 8 | 10718.9 | polar | flexible 4 | no |
| 9 | 21343.7 | azimuth | flexible 3 | no |
| 10 | 29518.9 | azimuth | flexible 4 | no |
| 11 | 39994.3 | polar | flexible 5 | no |
| 12 | 47898.7 | azimuth | flexible 5 | no |

States of the modal model: **14** (Guyan: 16). Modes added by the validation loop: **1**: azimuth mode 2 (4595.1 Hz).

### Comparison with the all-DOF model (T5 criteria, up to 2000 Hz)

Direct channels: largest deviation outside +/-2 % around the reference zeros; zeros: largest relative frequency deviation; 2x2: max ||G - G_ref||_2 / ||G_ref||_2.

| theta0 [deg] | Omega0 [deg/s] | modal: direct | modal: zeros | modal: 2x2 | modal T5 | Guyan: direct | Guyan: zeros | Guyan T5 |
|---|---|---|---|---|---|---|---|---|
| -45 | 0 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| -45 | 10 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| -15 | 0 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| -15 | 10 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 0 | 0 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 0 | 10 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 15 | 0 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 15 | 10 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 45 | 0 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |
| 45 | 10 | 0.0401 dB / 0.099 deg | 1.16e-04 | 2.20e-04 | yes | 22.1 dB / 157 deg | 1.44e-01 | no |

### Residual flexibility of the discarded modes and zeros of the isolated axis

| axis | collocated F [rad/(N m)] | non-collocated F [rad/(N m)] |
|---|---|---|
| azimuth | 7.144e-08 (9.4e-04 of the total) | 3.442e-09 (+4.5e-05) |
| polar | 6.799e-08 (3.4e-04 of the total) | 5.788e-10 (+6.1e-06) |

| axis | TF | f full [Hz] | pure truncation | truncation + static correction |
|---|---|---|---|---|
| azimuth | tau -> encoder | 12.3806 | 12.3809 (+2.25e-05) | 12.3806 (+4.03e-12) |
| azimuth | tau -> motor | 12.3806 | 12.3864 (+4.68e-04) | 12.3806 (+2.13e-11) |
| azimuth | tau -> motor | 1257.06 | 1258.03 (+7.78e-04) | 1257.06 (+3.66e-07) |
| polar | tau -> encoder | 16.9286 | 16.9287 (+2.83e-06) | 16.9286 (+7.10e-12) |
| polar | tau -> encoder | 504.488 | 504.429 (-1.16e-04) | 504.488 (-2.59e-07) |
| polar | tau -> motor | 11.8616 | 11.8636 (+1.68e-04) | 11.8616 (+1.06e-12) |
| polar | tau -> motor | 124.095 | 124.115 (+1.59e-04) | 124.095 (+1.10e-10) |
| polar | tau -> motor | 504.488 | 504.706 (+4.33e-04) | 504.488 (+4.95e-09) |

With `static_correction: residual_vectors` (residual static vectors at the motor and attachment nodes) the largest shift of the retained frequencies is azimuth: 3.0e-13, polar: 5.9e-10.

### Mesh check (lumped against consistent mass, retained modes)

| axis | mode | f consistent [Hz] | f lumped [Hz] | deviation | warning (> 2 %) |
|---|---|---|---|---|---|
| azimuth | 1 | 752.185 | 752.066 | -0.02 % | no |
| azimuth | 2 | 4595.13 | 4572.95 | -0.48 % | no |
| polar | 1 | 87.5063 | 87.5063 | -0.00 % | no |
| polar | 2 | 336.31 | 336.286 | -0.01 % | no |
| polar | 3 | 3506.35 | 3480.57 | -0.74 % | no |

### Main model (modal): operating points

| theta0 [deg] | Omega0 [deg/s] | tau_theta0 [N m] | first frequencies [Hz] |
|---|---|---|---|
| -45 | 0 | 7.2105 | 0.48378, 87.506, 336.31, 752.2, 3506.3, 4595.1 |
| -45 | 10 | 7.2013 | 0.48397, 87.506, 336.31, 752.2, 3506.3, 4595.1 |
| -15 | 0 | 10.401 | 0.31812, 87.506, 336.31, 752.19, 3506.3, 4595.1 |
| -15 | 10 | 10.397 | 0.31888, 87.506, 336.31, 752.19, 3506.3, 4595.1 |
| 0 | 0 | 10.977 | 0.14821, 87.506, 336.31, 752.18, 3506.3, 4595.1 |
| 0 | 10 | 10.978 | 0.14994, 87.506, 336.31, 752.18, 3506.3, 4595.1 |
| 15 | 0 | 10.805 | 87.506, 336.31, 752.19, 3506.3, 4595.1 |
| 15 | 10 | 10.81 | 87.506, 336.31, 752.19, 3506.3, 4595.1 |
| 45 | 0 | 8.3139 | 87.506, 336.31, 752.2, 3506.3, 4595.1 |
| 45 | 10 | 8.3231 | 87.506, 336.31, 752.2, 3506.3, 4595.1 |

## Nonlinear simulations

Model: modal (modal truncation, theta_ref from config.yaml).

- **S1** (doublet on tau_psi at Omega_0 = 0, theta0 = 0 deg, Omega0 = 0 deg/s): relative nonlinear/linear error psi_e: 3.22e-06, theta_e: 8.80e-01, psi_m: 3.85e-05, theta_m: 8.80e-01, psi_L: 6.05e-10, theta_p: 8.93e-01, theta_r: 8.79e-01, err_psi: 3.46e-06, err_theta: 3.49e-02. Peak of theta_e: 2.81e-06 arcsec.
- **S2** (tracking at Omega_max with a doublet on tau_psi, theta0 = 0 deg, Omega0 = 10 deg/s): relative nonlinear/linear error psi_e: 3.65e-06, theta_e: 8.82e-04, psi_m: 4.35e-05, theta_m: 8.82e-04, psi_L: 1.13e-08, theta_p: 8.82e-04, theta_r: 8.82e-04, err_psi: 3.91e-06, err_theta: 1.34e-03. Peak of theta_e: 0.00153 arcsec.
- **S3** (conservation with C = 0 and tau = 0): relative drift of E = 7.13e-15, of p_psi = 2.56e-12 (threshold 1e-8).
- **S1, comparison with the all-DOF model** (DOP853 with rtol 1e-8, atol 1e-10: the modes of the full FEM up to hundreds of kHz make the explicit integrator stability-limited; t_end = 0.1 s): maximum deviation relative to the peak psi_e: 1.64e-05, psi_m: 3.40e-04, psi_L: 5.97e-09.

In the polar channels of S1 the linear response exists only through the products of inertia: the quadratic term (1/2) J' dpsi^2 is of the same order there, so the relative error is not a measure of the quality of the linearisation.

## Lumped-parameter model (3 + 4 DOF)

> Every number in this section is a **result from placeholder parameters**.

Equivalent elements: J_eq = sum of the inertias, 1/k_eq = sum of the compliances (springs in series); masses lumped at the nodes (half of each adjacent element + node inertia). Coupling between the axes as in the coupled model (telescope at P, fork at L). Newton-Euler equations verified identical to the Lagrange ones.

| axis | element | nodes | FEM elements | J_eq [kg m^2] | k_eq [N m/rad] |
|---|---|---|---|---|---|
| azimuth | 1 | 1-5 | 1, 2, 3, 4 | 0.00036989 | 9547.55 |
| azimuth | 2 | 5-9 | 5, 6, 7, 8 | 5.4603e-05 | 13090.8 |
| polar | 1 | 1-H | 1, 2 | 1.7883e-05 | 37907.6 |
| polar | 2 | H-7 | 4, 5, 6 | 2.5094e-05 | 9630.88 |
| polar | 3 | 7-11 | 7, 8, 9, 10 | 0.0012883 | 9694.19 |
| polar | 4 | H-P | b | 0 | 9800 |

| axis | mode | lumped [Hz] | full FEM [Hz] |
|---|---|---|---|
| azimuth | 1 | 774.046 | 752.185 |
| azimuth | 2 | 1846.74 | 4595.13 |
| polar | 1 | 87.543 | 87.5063 |
| polar | 2 | 387.209 | 336.31 |
| polar | 3 | 992.253 | 3506.35 |
| polar | 4 | 10364.2 | 10718.9 |

Closed-form antiresonances against the zeros of the linearised coupled model (Hz):

- psi_e/tau_psi: 13.3162 (linear 13.3162)
- psi_m/tau_psi (1): 1143.53 (linear 1143.53)
- psi_m/tau_psi (2): 13.3162 (linear 13.3162)
- theta_e/tau_theta: 617.424 (linear 617.424)
- theta_e/tau_theta (tel): 16.9358 (linear 16.9358)

Against the all-DOF model (T5 criteria, theta0 = -45 deg): direct channels 25.4 dB / 160 deg, 2x2 matrix 18.7: the lumped model reproduces the rigid-body motion and the load-spring antiresonance, not the structural modes.

Simulations (nonlinear Newton-Euler simulator):

- S1 (theta0 = 0 deg, Omega0 = 0 deg/s): e(psi_e) = 1.32e-09, e(theta_p) = 8.93e-01
- S3: drift of E = 3.30e-13, of p_psi = 4.12e-15
- S4a (theta0 = -45 deg, Omega0 = 10 deg/s, offset 2 deg): e(psi_e) = 4.17e-03, e(theta_p) = 1.24e-02
- S4b (theta0 = -45 deg, Omega0 = 10 deg/s, offset 10 deg): e(psi_e) = 3.93e-02, e(theta_p) = 8.51e-02
- S4c (theta0 = 0 deg, Omega0 = 10 deg/s, offset 10 deg): e(psi_e) = 1.59e+01, e(theta_p) = 4.71e+00

Simulator code: `results_hub_telescopeonring_alhub_stiffring/simulatore_concentrato/ne_concentrato.py`, `results_hub_telescopeonring_alhub_stiffring/simulatore_concentrato/ne_concentrato.m`, `results_hub_telescopeonring_alhub_stiffring/simulatore_concentrato/ne_concentrato_parametri.m`.

## Comparison with the experimental measurements (7-DOF lumped model)

Kollmorgen WorkBench measurements (plant = encoder velocity / current, (rad/s)/A); Validation folder excluded.

| measurement | feature | measured f [Hz] | model f [Hz] | channel |
|---|---|---|---|---|
| azimuth, polar at 0 deg, 3.5 A | resonance | 3.797 | -- | -- |
| azimuth, polar at 0 deg, 3.5 A | antiresonance | 11.92 | 12.38 | psi_e |
| azimuth, polar at 0 deg, 3.5 A | resonance | 37.4 | -- | -- |
| azimuth, polar at 0 deg, 3.5 A | antiresonance | 45.26 | 87.54 | psi_e |
| azimuth, polar at 0 deg, 3.5 A | resonance | 490.4 | 774 | psi_e |
| azimuth, polar at 0 deg, 2.36 A | resonance | 4.177 | -- | -- |
| azimuth, polar at 0 deg, 2.36 A | antiresonance | 11.92 | 12.38 | psi_e |
| azimuth, polar at 0 deg, 2.36 A | resonance | 37.4 | -- | -- |
| azimuth, polar at 0 deg, 2.36 A | antiresonance | 45.26 | 87.54 | psi_e |
| azimuth, polar at 0 deg, 2.36 A | resonance | 490.4 | 774 | psi_e |
| azimuth, polar at 30 deg, 3.5 A | resonance | 3.82 | -- | -- |
| azimuth, polar at 30 deg, 3.5 A | antiresonance | 12.25 | 12.86 | psi_e |
| azimuth, polar at 30 deg, 3.5 A | resonance | 37.06 | -- | -- |
| azimuth, polar at 30 deg, 3.5 A | antiresonance | 44.14 | 87.54 | psi_e |
| azimuth, polar at 30 deg, 3.5 A | resonance | 481.3 | 774.1 | psi_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 8.4 | -- | -- |
| elevation, 3.5 A, 1-2000 Hz | resonance | 9.179 | -- | -- |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 18.12 | 16.93 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 152.2 | -- | -- |
| elevation, 3.5 A, 1-2000 Hz | resonance | 192.8 | -- | -- |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 300.3 | -- | -- |
| elevation, 3.5 A, 1-2000 Hz | resonance | 348.2 | 387.2 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 358.6 | 617.4 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 380.4 | 387.2 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 403.6 | 617.4 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 428.2 | 387.2 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 454.3 | 387.2 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 526.6 | 617.4 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 667.1 | 992.3 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 707.7 | 992.3 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | resonance | 773.3 | 992.3 | theta_e |
| elevation, 3.5 A, 1-2000 Hz | antiresonance | 1039 | 774 | theta_e |
| elevation without telescope, 2.5 A | antiresonance | 404.1 | 617.4 | theta_e |
| elevation without telescope, 2.5 A | resonance | 428.4 | 393.1 | theta_e |
| elevation without telescope, 2.5 A | antiresonance | 510.2 | 617.4 | theta_e |

- Azimuth level 15-300 Hz: measurement -4.9 dB relative to the model.
- Equivalent stiffness to the ground (azimuth): k_0 ~ 1.23e+03 N m/rad (absent from the model).
- Motor-side inertia for the first measured mode: 0.00099 kg m^2 (model 0.000397).
- Without telescope: k_0 ~ 2.75e+03 N m/rad, J ~ 0.0358 kg m^2 (model 0.0183).
- Elevation: antiresonance at the encoder 18.12 Hz, resonance 48 Hz; two-mass fit J_1 ~ 0.13 kg m^2, J_c ~ 0.77 kg m^2, k_b ~ 1e+04 N m/rad (telescope hanging from the shaft through a compliance; the model has no zero of theta_e below 16.93 Hz).

## Figures

- `results_hub_telescopeonring_alhub_stiffring/figures/geometria_azimut.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/modi_azimut.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/confronto_completo_ridotto_azimut.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/geometria_polare.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/modi_polare.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/confronto_completo_ridotto_polare.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/geometria_mozzo.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th-45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th-45_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/errore_misura_th-45_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th-15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th-15_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/errore_misura_th-15_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+0_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+0_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/errore_misura_th+0_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+15_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/errore_misura_th+15_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_encoder_th+45_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/errore_misura_th+45_om10.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_modale_th-45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_modale_th-15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_modale_th+0_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_modale_th+15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_modale_th+45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/modi_modale_azimut.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/modi_modale_polare.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/frequenze_vs_theta0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/frf_incrociate_vs_omega.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_S1.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_S2.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_S3.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_concentrato_th-45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_concentrato_th-15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_concentrato_th+0_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_concentrato_th+15_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/bode_concentrato_th+45_om0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/frequenze_vs_theta0_concentrato.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/catene_polari_poli_zeri.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_concentrato_S4a.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_concentrato_S4b.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_concentrato_S4c.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/simulazione_concentrato_S1.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_Azimuth_BodeMeasurement_Po_0.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_Azimuth_BodeMeasurement_Po_0_dep.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_Azimuth_BodeMeasurement_Po_30.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_Elevation_BodeMeasurement_bode_New_16_10_26.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_BodeMeasurement_NoTelescope.pdf`
- `results_hub_telescopeonring_alhub_stiffring/figures/exp_azimut_sovrapposte.pdf`

## Tables

- `results_hub_telescopeonring_alhub_stiffring/tables/element_breakdown_azimut.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/element_breakdown_polare.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/modal_modes_azimut.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/modal_modes_polare.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/modal_spectrum.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/modal_vs_full.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/reduction_validation.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/full_model_frequencies.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/poles_zeros.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/operating_points.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/lumped_parameters.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/lumped_vs_full.csv`
- `results_hub_telescopeonring_alhub_stiffring/tables/experimental_comparison.csv`

## LaTeX document

`results_hub_telescopeonring_alhub_stiffring/latex/multidof_fem_model.tex`: compilato, PDF `results_hub_telescopeonring_alhub_stiffring/latex/multidof_fem_model.pdf`

## Status of the decisions and open questions (section 3)

| ID | status |
|---|---|
| D1 | Confirmed: telescope on polar element 3, node P at the midpoint. |
| D2 | Confirmed: polar elements 2 and 4 on the bearings. |
| D3 | Confirmed: polar masters [1, P, 7, 11]. |
| Q1 | **Default: xi_p = 0.5 (telescope over the whole element 3).** (assumption) |
| Q2 | **Default: no stiffening by the clamp (payload_clamp_rigid_mm = 0).** (assumption) |
| Q3 | User decision (2026-10-12): telescope tensor entered in config.yaml (frame `cad_cg`), subtraction of the polar parts: no. |
| Q4 | User decision (2026-10-13): Gawronski modal truncation of the full FEM (reduction.method: modal); Guyan (decision of 2026-10-12) and Craig-Bampton remain as comparisons. |
| Q5 | **Default: stiffnesses from the geometry, no joint compliance.** (assumption) |
| Q6 | User answer (2026-10-12): values of the .m file (steel 8000 kg/m^3, 77 GPa; aluminium 7075 2810 kg/m^3, 26.9 GPa). Still placeholders per part. |
| Q7 | **Default: null, J^ax without Steiner terms (irrelevant if the azimuth load includes the polar shafts).** (assumption) |
| Q8 | **Default: C = beta K with zeta = 0.01 on the first flexible mode; bearings b = 0.** (assumption) |
| Q9 | User data (2026-10-13): load 'load (fork + polar shafts)' from Onshape (Onshape, Tavola Rotante [DEF] <1>, without mate connector), m = 27.831 kg, J = I_zz,G + m d^2 = 1.2480 + 27.831 x 12.212^2 mm^2 = 1.25215 kg m^2; it includes the polar shafts: J^ax not added (already in the CAD value). |
| Q10 | Solved: z_O = 54.9 cm measured in the CAD and already confirmed in the 2-DOF model. |
| Q11 | **Default: Omega_max = 10 deg/s.** (assumption) |
| Q12 | Assumption adopted and declared: rigid fork. |

## Warnings

- **[PLACEHOLDER] Q6_MATERIALI**: Q6: placeholder materials ['steel', 'aluminium'] (steel: rho = 8000 kg/m^3, G = 77 GPa, aluminium: rho = 2810 kg/m^3, G = 26.9 GPa).
- **[NOTE] TELESCOPIO_SENZA_ALBERO**: Polar element 3 in mode: telescope: no shaft under the telescope. In the hub model its nodes 3 and 4 coincide with the hub centre H; the telescope is the node P attached to H by k_b, and its inertia enters only through the payload section.
- **[PLACEHOLDER] MOZZO_TRATTO_TELESCOPIO**: hub.telescope_on_ring: the section of the telescope tube clamped in the support ring (300/290 x 100 mm, alluminio) is assumed rigid with the ring: its inertia is moved from the telescope to the hub centre H.
- **[NOTE] MOZZO_TELESCOPIO**: Hub model (polar.hub.enabled: true): the hub (two side hubs, tube + flange, and the support ring) replaces polar elements [2, 4] and its centre H takes the nodes of element 3; the telescope is a node P off the chain, attached to H by k_b = 9800 N m/rad.
- **[PLACEHOLDER] MOZZO_ANELLO_IRRIGIDITO**: hub.ring.stiffness_factor = 100: the stiffness of the support ring from the half-ring model (ring alone, lower bound) is multiplied by this factor, chosen from the test without the telescope (no structural mode below ~400 Hz).
- **[PLACEHOLDER] MOZZO_K_B**: k_b = 9800 N m/rad (telescope to hub) is a placeholder (placeholder: two-mass estimate from the antiresonance (18.1 Hz) and resonance (~48 Hz) measured at the elevation encoder): it sets the antiresonance at the encoder.
- **[PLACEHOLDER] MOZZO_INERZIA_AZIMUT**: hub.inertia.in_azimuth_load = false: the inertia of the hub about z_hat is added to J^ax. Check whether the azimuth load exported from Onshape already contains the hub parts.
- **[PLACEHOLDER] NODI_SENZA_INERZIA_AZIMUT**: Azimuth axis: no concentrated inertia at encoder (node 1), motor (node 5) (rotor/disc missing): those nodes carry only the shaft mass and are almost inertia-free. This affects the zeros and the high modes (typically strongly) and the choice of the truncation modes.
- **[PLACEHOLDER] Q5_GIUNZIONI_AZIMUT**: Q5: azimuth axis without joint compliances (shrink fits, keys, flanges): the stiffnesses are those of the shafts alone and probably overestimate the real one.
- **[PLACEHOLDER] Q8_CUSCINETTI_AZIMUT**: Q8: zero viscous bearing friction on the azimuth axis.
- **[PLACEHOLDER] NODI_SENZA_INERZIA_POLARE**: Polar axis: no concentrated inertia at encoder (node 1), motor (node 7) (node_inertias entirely commented out): those nodes carry only the shaft mass and are almost inertia-free. This affects the zeros and the high modes (typically strongly) and the choice of the truncation modes.
- **[PLACEHOLDER] Q5_GIUNZIONI_POLARE**: Q5: polar axis without joint compliances (shrink fits, keys, flanges): the stiffnesses are those of the shafts alone and probably overestimate the real one.
- **[NOTE] Q7_POSIZIONE_POLARE**: Q7: y_node1_from_o2_mm = null. The contribution of the polar shaft to the azimuth inertia (J^ax) is computed from the transverse moments of the parts about their own centres only, without Steiner terms about o_2. Irrelevant: the azimuth load includes the polar shafts (includes_polar_shafts: true), so J^ax is not added.
- **[PLACEHOLDER] Q8_CUSCINETTI_POLARE**: Q8: zero viscous bearing friction on the polar axis.
- **[NOTE] TERNE_ORIGINI**: Origins on z_hat: dist(O_N, o_1) = 59.695 mm, dist(O_N, o_2) = 548.865 mm (o_2 = centre of rotation of the polar shaft). O_N is the origin of the Onshape assembly (centre of the base of the azimuth shaft) and x_hat, y_hat, z_hat are its axes.
- **[NOTE] RIDUZIONE_MODALE**: Reduction: Gawronski modal truncation of the full FEM of each axis (band 2000 Hz, cut-off 4000 Hz, static correction: none).
- **[PLACEHOLDER] Q8_SMORZAMENTO**: Q8: stiffness-proportional damping, C = beta K, with zeta = 0.01 on the first flexible mode of each axis. Placeholder: the height of the FRF peaks and notches is not representative.
- **[PLACEHOLDER] Q11_OMEGA_MAX**: Q11: tracking rates [0.0, 10.0] deg/s; Omega_max = 10 deg/s placeholder. It sets the gyroscopic coupling.
- **[NOTE] Q12_FORCELLA_RIGIDA**: Q12: rigid fork between azimuth node 9 and the polar bearings. The reaction of tau_theta does no work on psi_L (section 9.4).
- **[PLACEHOLDER] MISURA_THETA0_3**: Elevation_BodeMeasurement_bode_New_16_10_26.csv: elevation during the test not recorded in the file; theta0 = 0 deg assumed.
- **[NOTE] MISURE_SPERIMENTALI**: Comparison with 5 plant measurements (Kollmorgen WorkBench, LucyIdentification); K_t = 0.258 N m/A. The Validation folder is not used.
- **[NOTE] JAX_NEL_CARICO**: The azimuth load 'load (fork + polar shafts)' (node 9) includes the polar shafts (includes_polar_shafts: true): J^ax = 8.3476e-04 kg m^2 computed from the polar elements is not added, because their moment about the azimuth axis, with the real distances from the axis, is already in the CAD value.
- **[PLACEHOLDER] PAYLOAD_SEGNO_PRODOTTI**: payload.products_sign = 'tensor': the sign convention of the Onshape products of inertia is not verified on a test body; it affects the products I_IJ, I_IK, I_JK, not the diagonal moments.
- **[NOTE] PAYLOAD_TRATTO_SU_MOZZO**: Telescope tensor without the tube section clamped in the support ring: m - 1.302 kg, I about the pointing axis - 2.8337e-02, about the polar axis and z - 1.5254e-02 kg m^2 (moved to the hub centre H).
- **[NOTE] Q3_TELESCOPIO**: Q3: telescope from the tensor entered in config.yaml (frame: cad_cg, Onshape, Telescopio LAMERA <1>, without mate connector), without subtractions (subtract_shaft_parts: false): the tensor must contain the telescope only, not the parts of the polar axis already modelled in the elements.

