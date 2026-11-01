# Grey-box estimate of H7 on the VTesi plant responses (config_vtesi_greybox.yaml)

Model `../config_hub_telescopeonring_alhub_fit_ground_aztel.yaml`, estimator `lucy_cl.greyest`: error log, loss soft_l1 (f_scale 0.5), band 1-400 Hz, covariance robust. Bounds: inertias and geometric springs +-20 %, the other springs and the damping free.

## Azimuth

Measurement `BodeMeasurement_Azimuth_14viti.csv`; best start {'k_{\\psi,3}': 300000.0}; converged, 48 evaluations; loss 293.3 -> 103.4; Jacobian condition 6.65e+06.

| fit [%] | all | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| H7 nominal | 35.7 | -72.6 | -11.4 | 18.8 |
| grey-box estimate | 51.2 | 75.4 | 22.3 | 34.7 |

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\psi,e}` | inertia | 0.0001849 | 0.000148 - 0.000222 | 0.0002219 **(at bound)** | 1.2 | 1899.8% |  |
| `J_{\psi,m}` | inertia | 0.0002941 | 0.000235 - 0.000353 | 0.0002353 **(at bound)** | 0.8 | 2178.1% | rotor included |
| `J_{\psi,L}` | inertia | 1.252 | 1 - 1.5 | 1.466 | 1.17 | 28.7% |  |
| `J^{ax}` | inertia | 0.05331 | 0.0426 - 0.064 | 0.06397 **(at bound)** | 1.2 | 1.3% | hub and ring about z |
| `I_zz` | inertia | 0.8578 | 0.686 - 1.03 | 1.029 **(at bound)** | 1.2 | 40.3% | telescope about z at theta = 0 |
| `k_{\psi,1}` | geometric | 9548 | 7.64e+03 - 1.15e+04 | 7638 **(at bound)** | 0.8 | 1826.9% | shafts |
| `k_{\psi,2}` | geometric | 1.309e+04 | 1.05e+04 - 1.57e+04 | 1.571e+04 **(at bound)** | 1.2 | 2.9% | shafts |
| `k_{\psi,3}` | free | 1.16e+05 | 1 - 1e+09 | 8.105e+04 | 0.699 | 12.4% | telescope - fork |
| `c_{\psi,1}` | damping | 0.4034 | 1e-06 - 1e+04 | 1e-06 **(at bound)** | 2.48e-06 | 0.0% |  |
| `c_{\psi,2}` | damping | 0.5531 | 1e-06 - 1e+04 | 5.111 | 9.24 | 10.4% |  |
| `c_{\psi,3}` | damping | 4.901 | 1e-06 - 1e+04 | 5.801 | 1.18 | 109.7% |  |
| `kg_az` | free | 1418 | 1 - 1e+09 | 1043 | 0.735 | 1.9% | ground constraint stiffness |
| `eta_az` | damping | 0.32 | 0.0001 - 100 | 0.7255 | 2.27 | 2.5% | ground constraint loss factor |

Starts: {'k_{\\psi,3}': 30000.0} -> loss 103.4; {'k_{\\psi,3}': 60000.0} -> loss 103.4; {'k_{\\psi,3}': 120000.0} -> loss 103.4; {'k_{\\psi,3}': 300000.0} -> loss 103.4

Not determined by the data (std above 50 % of the estimate): `J_{\psi,e}`, `J_{\psi,m}`, `k_{\psi,1}`, `c_{\psi,3}`.

## Elevation

Measurement `BodeMeasurement_Polar_14viti.csv`; best start {'k_{\\theta,4}': 400000.0, 'ring': 1.0}; converged, 95 evaluations; loss 1625 -> 704.7; Jacobian condition 8.18e+07.

| fit [%] | all | 1-5 Hz | 5-100 Hz | >100 Hz |
|---|---|---|---|---|
| H7 nominal | -511.4 | 48.9 | 19.4 | -5409.9 |
| grey-box estimate | 43.3 | 41.5 | 18.0 | -8.4 |

| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |
|---|---|---|---|---|---|---|---|
| `J_{\theta,e}` | inertia | 8.942e-06 | 7.15e-06 - 1.07e-05 | 1.073e-05 **(at bound)** | 1.2 | 505565.2% |  |
| `J_{\theta,h}` | inertia | 0.03227 | 0.0258 - 0.0387 | 0.02582 **(at bound)** | 0.8 | 29220.7% |  |
| `J_{\theta,m}` | inertia | 0.0007386 | 0.000591 - 0.000886 | 0.0005911 | 0.8 | 472619.4% | rotor included |
| `J_{\theta,r}` | inertia | 0.0006441 | 0.000515 - 0.000773 | 0.0005153 **(at bound)** | 0.8 | 526729.9% |  |
| `I_yy` | inertia | 0.8662 | 0.693 - 1.04 | 0.693 **(at bound)** | 0.8 | 296.3% | telescope about the polar axis |
| `k_{\theta,1}^{geo}` | geometric | 4.798e+04 | 3.84e+04 - 5.76e+04 | 3.838e+04 **(at bound)** | 0.8 | 55603.2% | hub and shafts, ring excluded |
| `k_{\theta,2}^{geo}` | geometric | 1.017e+04 | 8.14e+03 - 1.22e+04 | 1.22e+04 | 1.2 | 381289.6% | hub and shafts, ring excluded |
| `k_{\theta,3}` | geometric | 9694 | 7.76e+03 - 1.16e+04 | 1.163e+04 **(at bound)** | 1.2 | 563987.3% | shafts |
| `k_{\theta,4}` | free | 1.16e+05 | 1 - 1e+09 | 15.57 | 0.000134 | 1260.1% | telescope - hub |
| `ring` | free | 1000 | 0.001 - 1e+06 | 2.625 | 0.00262 | 531031.4% | stiffness factor of the support ring |
| `c_{\theta,1}` | damping | 0.5053 | 1e-06 - 1e+04 | 3.887e-06 | 7.69e-06 | 63.3% |  |
| `c_{\theta,2}` | damping | 0.1094 | 1e-06 - 1e+04 | 23.13 | 211 | 488486.1% |  |
| `c_{\theta,3}` | damping | 0.1048 | 1e-06 - 1e+04 | 3.192e-06 | 3.04e-05 | 284.0% |  |
| `c_{\theta,4}` | damping | 1.254 | 1e-06 - 1e+04 | 483.7 | 386 | 606.7% |  |
| `kg_el` | free | 164 | 1 - 1e+09 | 106.7 | 0.651 | 14.2% | ground constraint stiffness |
| `eta_el` | damping | 3.11 | 0.0001 - 100 | 4.638 | 1.49 | 14.5% | ground constraint loss factor |

Starts: {'k_{\\theta,4}': 60000.0, 'ring': 1.0} -> loss 715.5; {'k_{\\theta,4}': 60000.0, 'ring': 1000.0} -> loss 720.5; {'k_{\\theta,4}': 150000.0, 'ring': 1.0} -> loss 715.5; {'k_{\\theta,4}': 150000.0, 'ring': 1000.0} -> loss 715.5; {'k_{\\theta,4}': 400000.0, 'ring': 1.0} -> loss 704.7; {'k_{\\theta,4}': 400000.0, 'ring': 1000.0} -> loss 704.9; {'k_{\\theta,4}': 1000000.0, 'ring': 1.0} -> loss 910.8; {'k_{\\theta,4}': 1000000.0, 'ring': 1000.0} -> loss 963.2

Not determined by the data (std above 50 % of the estimate): `J_{\theta,e}`, `J_{\theta,h}`, `J_{\theta,m}`, `J_{\theta,r}`, `I_yy`, `k_{\theta,1}^{geo}`, `k_{\theta,2}^{geo}`, `k_{\theta,3}`, `k_{\theta,4}`, `ring`, `c_{\theta,1}`, `c_{\theta,2}`, `c_{\theta,3}`, `c_{\theta,4}`.

