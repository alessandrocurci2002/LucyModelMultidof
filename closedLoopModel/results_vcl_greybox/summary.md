# Grey-box estimate: H7 (9 DOF) estimated on the VClosedLoop plant responses (20 Oct 2026) (config_vcl_greybox.yaml)

Model `../config_hub_telescopeonring_alhub_fit_ground_aztel.yaml`, estimator `lucy_cl.greyest`: error log, loss soft_l1 (f_scale 0.5), band 1-150 Hz, covariance robust. Bounds: inertias and geometric springs +-20 %, the other springs free; damping limited to zeta <= 0.1 on every mode up to inf Hz (ground loss factor excluded).

## Azimuth

Measurement `BodeMeasurement_azimuth.csv`; best start {'k_{\\psi,3}': 60000.0}; converged, 68 evaluations; loss 1.424e+04 -> 13.01; Jacobian condition 1.98e+07.

| fit [%] | band 1-150 Hz | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| nominal | -30.7 | -126.5 | -7.4 | -94.8 |
| grey-box estimate | 53.8 | 84.9 | 22.6 | -100.5 |

Modes of the axis in the estimated model up to inf Hz (ground stiffness included; zeta 1 = overdamped): 1.4 Hz (zeta 0.000), 55.0 Hz (zeta 0.004), 899.9 Hz (zeta 0.100), 1991.4 Hz (zeta 0.079).

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\psi,e}` | inertia | 0.0001849 | 0.000148 - 0.000222 | 0.000148 **(at bound)** | 0.8 | 182713.9% |  |
| `J_{\psi,m}` | inertia | 0.0002941 | 0.000235 - 0.000353 | 0.0002353 **(at bound)** | 0.8 | 64837.0% | rotor included |
| `J_{\psi,L}` | inertia | 1.252 | 1 - 1.5 | 1.495 | 1.19 | 82.3% |  |
| `J^{ax}` | inertia | 0.05331 | 0.0426 - 0.064 | 0.06366 | 1.19 | 657.8% | hub and ring about z |
| `I_zz` | inertia | 0.8578 | 0.686 - 1.03 | 1.029 | 1.2 | 104.1% | telescope about z at theta = 0 |
| `k_{\psi,1}` | geometric | 9548 | 7.64e+03 - 1.15e+04 | 1.146e+04 **(at bound)** | 1.2 | 244009.3% | shafts |
| `k_{\psi,2}` | geometric | 1.309e+04 | 1.05e+04 - 1.57e+04 | 1.521e+04 | 1.16 | 11.6% | shafts |
| `k_{\psi,3}` | free | 1.16e+05 | 1 - 1e+09 | 7.398e+04 | 0.638 | 33.5% | telescope - fork |
| `c_{\psi,1}` | damping | 0.4034 | 1e-06 - 1e+04 | 1.318e-06 | 3.27e-06 | 0.8% |  |
| `c_{\psi,2}` | damping | 0.5531 | 1e-06 - 1e+04 | 0.729 | 1.32 | 354.5% |  |
| `c_{\psi,3}` | damping | 4.901 | 1e-06 - 1e+04 | 1.527 | 0.311 | 184.7% |  |
| `kg_az` | free | 1418 | 1 - 1e+09 | 214.7 | 0.151 | 5.7% | ground constraint stiffness |
| `eta_az` | damping | 0.32 | 0.0001 - 100 | 1.221 | 3.81 | 6.1% | ground constraint loss factor |

Starts: {'k_{\\psi,3}': 30000.0} -> loss 14.78; {'k_{\\psi,3}': 60000.0} -> loss 13.01; {'k_{\\psi,3}': 120000.0} -> loss 20.62; {'k_{\\psi,3}': 300000.0} -> loss 21.84

Not determined by the data (std above 50 % of the estimate): `J_{\psi,e}`, `J_{\psi,m}`, `J_{\psi,L}`, `J^{ax}`, `I_zz`, `k_{\psi,1}`, `c_{\psi,2}`, `c_{\psi,3}`.

## Elevation

Measurement `BodeMeasurement_elevation.csv`; best start {'k_{\\theta,4}': 400000.0, 'ring': 1000.0}; converged, 74 evaluations; loss 8.483e+04 -> 110.5; Jacobian condition 1.11e+08.

| fit [%] | band 1-150 Hz | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| nominal | 46.3 | 51.1 | 25.0 | -922.6 |
| grey-box estimate | 44.4 | 42.6 | 22.1 | 4.4 |

Modes of the axis in the estimated model up to inf Hz (ground stiffness included; zeta 1 = overdamped): 2.1 Hz (zeta 0.000), 155.7 Hz (zeta 0.052), 710.8 Hz (zeta 0.100), 2165.5 Hz (zeta 0.000), 3075.4 Hz (zeta 0.027).

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\theta,e}` | inertia | 8.942e-06 | 7.15e-06 - 1.07e-05 | 1.073e-05 | 1.2 | 35898689.0% |  |
| `J_{\theta,h}` | inertia | 0.03227 | 0.0258 - 0.0387 | 0.02999 | 0.929 | 70371780.4% |  |
| `J_{\theta,m}` | inertia | 0.0007386 | 0.000591 - 0.000886 | 0.0008842 | 1.2 | 8480998.8% | rotor included |
| `J_{\theta,r}` | inertia | 0.0006441 | 0.000515 - 0.000773 | 0.0007699 | 1.2 | 9093737.1% |  |
| `I_yy` | inertia | 0.8662 | 0.693 - 1.04 | 0.7201 | 0.831 | 2930720.8% | telescope about the polar axis |
| `k_{\theta,1}^{geo}` | geometric | 4.798e+04 | 3.84e+04 - 5.76e+04 | 4.526e+04 | 0.943 | 1582275.9% | hub and shafts, ring excluded |
| `k_{\theta,2}^{geo}` | geometric | 1.017e+04 | 8.14e+03 - 1.22e+04 | 8146 | 0.801 | 10623888.8% | hub and shafts, ring excluded |
| `k_{\theta,3}` | geometric | 9694 | 7.76e+03 - 1.16e+04 | 7833 | 0.808 | 13426476.2% | shafts |
| `k_{\theta,4}` | free | 1.16e+05 | 1 - 1e+09 | 1.074e+07 | 92.6 | 5835751.7% | telescope - hub |
| `ring` | free | 1000 | 0.001 - 1e+06 | 1.151 | 0.00115 | 2822098.3% | stiffness factor of the support ring |
| `c_{\theta,1}` | damping | 0.5053 | 1e-06 - 1e+04 | 1.069e-06 | 2.11e-06 | 813.7% |  |
| `c_{\theta,2}` | damping | 0.1094 | 1e-06 - 1e+04 | 0.1821 | 1.67 | 189473.6% |  |
| `c_{\theta,3}` | damping | 0.1048 | 1e-06 - 1e+04 | 0.3248 | 3.1 | 13625169.6% |  |
| `c_{\theta,4}` | damping | 1.254 | 1e-06 - 1e+04 | 30.39 | 24.2 | 56659203.6% |  |
| `kg_el` | free | 164 | 1 - 1e+09 | 128.9 | 0.786 | 35.6% | ground constraint stiffness |
| `eta_el` | damping | 3.11 | 0.0001 - 100 | 4.139 | 1.33 | 36.4% | ground constraint loss factor |

Starts: {'k_{\\theta,4}': 60000.0, 'ring': 1.0} -> loss 110.6; {'k_{\\theta,4}': 60000.0, 'ring': 1000.0} -> loss 139.2; {'k_{\\theta,4}': 150000.0, 'ring': 1.0} -> loss 139.1; {'k_{\\theta,4}': 150000.0, 'ring': 1000.0} -> loss 139.1; {'k_{\\theta,4}': 400000.0, 'ring': 1.0} -> loss 149; {'k_{\\theta,4}': 400000.0, 'ring': 1000.0} -> loss 110.5; {'k_{\\theta,4}': 1000000.0, 'ring': 1.0} -> loss 196.8; {'k_{\\theta,4}': 1000000.0, 'ring': 1000.0} -> loss 110.6

Not determined by the data (std above 50 % of the estimate): `J_{\theta,e}`, `J_{\theta,h}`, `J_{\theta,m}`, `J_{\theta,r}`, `I_yy`, `k_{\theta,1}^{geo}`, `k_{\theta,2}^{geo}`, `k_{\theta,3}`, `k_{\theta,4}`, `ring`, `c_{\theta,1}`, `c_{\theta,2}`, `c_{\theta,3}`, `c_{\theta,4}`.

