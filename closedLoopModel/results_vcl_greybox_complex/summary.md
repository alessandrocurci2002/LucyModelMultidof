# Grey-box estimate: H7 (9 DOF) estimated on the VClosedLoop plant responses (20 Oct 2026), complex error (config_vcl_greybox_complex.yaml)

Model `../config_hub_telescopeonring_alhub_fit_ground_aztel.yaml`, estimator `lucy_cl.greyest`: error complex, loss linear (f_scale 1.0), band 1-150 Hz, covariance robust. Bounds: inertias and geometric springs +-20 %, the other springs free; damping limited to zeta <= 0.1 on every mode up to inf Hz (ground loss factor excluded).

## Azimuth

Measurement `BodeMeasurement_azimuth.csv`; best start {'k_{\\psi,3}': 60000.0}; maximum evaluations reached, 3000 evaluations; loss 1.409e+04 -> 0.0004162; Jacobian condition 2.69e+08.

| fit [%] | band 1-150 Hz | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| nominal | -30.7 | -126.5 | -7.4 | -94.8 |
| grey-box estimate | 57.7 | 73.0 | 38.0 | -151.1 |

Modes of the axis in the estimated model up to inf Hz (ground stiffness included; zeta 1 = overdamped): 1.3 Hz (zeta 0.000), 51.9 Hz (zeta 0.016), 772.5 Hz (zeta 0.045), 1679.5 Hz (zeta 0.100).

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\psi,e}` | inertia | 0.0001849 | 0.000148 - 0.000222 | 0.0001866 | 1.01 | 608.2% |  |
| `J_{\psi,m}` | inertia | 0.0002941 | 0.000235 - 0.000353 | 0.0002974 | 1.01 | 344.7% | rotor included |
| `J_{\psi,L}` | inertia | 1.252 | 1 - 1.5 | 1.002 | 0.8 | 26.4% |  |
| `J^{ax}` | inertia | 0.05331 | 0.0426 - 0.064 | 0.04265 **(at bound)** | 0.8 | 1.1% | hub and ring about z |
| `I_zz` | inertia | 0.8578 | 0.686 - 1.03 | 1.029 | 1.2 | 26.1% | telescope about z at theta = 0 |
| `k_{\psi,1}` | geometric | 9548 | 7.64e+03 - 1.15e+04 | 9660 | 1.01 | 450.8% | shafts |
| `k_{\psi,2}` | geometric | 1.309e+04 | 1.05e+04 - 1.57e+04 | 1.506e+04 | 1.15 | 12.8% | shafts |
| `k_{\psi,3}` | free | 1.16e+05 | 1 - 1e+09 | 5.519e+04 | 0.476 | 4.7% | telescope - fork |
| `c_{\psi,1}` | damping | 0.4034 | 1e-06 - 1e+04 | 0.1872 | 0.464 | 190.3% |  |
| `c_{\psi,2}` | damping | 0.5531 | 1e-06 - 1e+04 | 0.2713 | 0.49 | 518.7% |  |
| `c_{\psi,3}` | damping | 4.901 | 1e-06 - 1e+04 | 5.377 | 1.1 | 61.5% |  |
| `kg_az` | free | 1418 | 1 - 1e+09 | 147 | 0.104 | 9.4% | ground constraint stiffness |
| `eta_az` | damping | 0.32 | 0.0001 - 100 | 1.833 | 5.73 | 9.5% | ground constraint loss factor |

Starts: {'k_{\\psi,3}': 30000.0} -> loss 0.0004635; {'k_{\\psi,3}': 60000.0} -> loss 0.0004162; {'k_{\\psi,3}': 120000.0} -> loss 0.000807; {'k_{\\psi,3}': 300000.0} -> loss 0.0005995

Not determined by the data (std above 50 % of the estimate): `J_{\psi,e}`, `J_{\psi,m}`, `k_{\psi,1}`, `c_{\psi,1}`, `c_{\psi,2}`, `c_{\psi,3}`.

## Elevation

Measurement `BodeMeasurement_elevation.csv`; best start {'k_{\\theta,4}': 150000.0, 'ring': 1000.0}; converged, 251 evaluations; loss 8.462e+04 -> 0.000334; Jacobian condition 2.3e+06.

| fit [%] | band 1-150 Hz | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| nominal | 46.3 | 51.1 | 25.0 | -922.6 |
| grey-box estimate | 54.5 | 60.4 | 27.8 | -1.0 |

Modes of the axis in the estimated model up to inf Hz (ground stiffness included; zeta 1 = overdamped): 2.1 Hz (zeta 0.000), 311.1 Hz (zeta 0.060), 781.2 Hz (zeta 0.096), 11577.8 Hz (zeta 0.100), 13007.9 Hz (zeta 0.000).

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\theta,e}` | inertia | 8.942e-06 | 7.15e-06 - 1.07e-05 | 9.109e-06 | 1.02 | 2439598.6% |  |
| `J_{\theta,h}` | inertia | 0.03227 | 0.0258 - 0.0387 | 0.03294 | 1.02 | 36527917.0% |  |
| `J_{\theta,m}` | inertia | 0.0007386 | 0.000591 - 0.000886 | 0.0008863 **(at bound)** | 1.2 | 4407728.6% | rotor included |
| `J_{\theta,r}` | inertia | 0.0006441 | 0.000515 - 0.000773 | 0.000773 **(at bound)** | 1.2 | 2750966.3% |  |
| `I_yy` | inertia | 0.8662 | 0.693 - 1.04 | 0.9877 | 1.14 | 1218368.0% | telescope about the polar axis |
| `k_{\theta,1}^{geo}` | geometric | 4.798e+04 | 3.84e+04 - 5.76e+04 | 4.888e+04 | 1.02 | 3493401.3% | hub and shafts, ring excluded |
| `k_{\theta,2}^{geo}` | geometric | 1.017e+04 | 8.14e+03 - 1.22e+04 | 8139 **(at bound)** | 0.8 | 1433691.7% | hub and shafts, ring excluded |
| `k_{\theta,3}` | geometric | 9694 | 7.76e+03 - 1.16e+04 | 7755 **(at bound)** | 0.8 | 5744529.8% | shafts |
| `k_{\theta,4}` | free | 1.16e+05 | 1 - 1e+09 | 2.128e+08 | 1.83e+03 | 4913208.9% | telescope - hub |
| `ring` | free | 1000 | 0.001 - 1e+06 | 2010 | 2.01 | 8334652.9% | stiffness factor of the support ring |
| `c_{\theta,1}` | damping | 0.5053 | 1e-06 - 1e+04 | 0.1325 | 0.262 | 34938345.0% |  |
| `c_{\theta,2}` | damping | 0.1094 | 1e-06 - 1e+04 | 0.5982 | 5.47 | 1575631.7% |  |
| `c_{\theta,3}` | damping | 0.1048 | 1e-06 - 1e+04 | 0.2088 | 1.99 | 11946562.7% |  |
| `c_{\theta,4}` | damping | 1.254 | 1e-06 - 1e+04 | 0.2847 | 0.227 | 28367769.3% |  |
| `kg_el` | free | 164 | 1 - 1e+09 | 173.1 | 1.06 | 12.2% | ground constraint stiffness |
| `eta_el` | damping | 3.11 | 0.0001 - 100 | 2.498 | 0.803 | 11.9% | ground constraint loss factor |

Starts: {'k_{\\theta,4}': 60000.0, 'ring': 1.0} -> loss 1.62e+06; {'k_{\\theta,4}': 60000.0, 'ring': 1000.0} -> loss 0.0003343; {'k_{\\theta,4}': 150000.0, 'ring': 1.0} -> loss 1.62e+06; {'k_{\\theta,4}': 150000.0, 'ring': 1000.0} -> loss 0.000334; {'k_{\\theta,4}': 400000.0, 'ring': 1.0} -> loss 1.62e+06; {'k_{\\theta,4}': 400000.0, 'ring': 1000.0} -> loss 0.0003355; {'k_{\\theta,4}': 1000000.0, 'ring': 1.0} -> loss 1.62e+06; {'k_{\\theta,4}': 1000000.0, 'ring': 1000.0} -> loss 0.0003347

Not determined by the data (std above 50 % of the estimate): `J_{\theta,e}`, `J_{\theta,h}`, `J_{\theta,m}`, `J_{\theta,r}`, `I_yy`, `k_{\theta,1}^{geo}`, `k_{\theta,2}^{geo}`, `k_{\theta,3}`, `k_{\theta,4}`, `ring`, `c_{\theta,1}`, `c_{\theta,2}`, `c_{\theta,3}`, `c_{\theta,4}`.

