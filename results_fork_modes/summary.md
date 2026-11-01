# Masse modali sullo statore polare e load di azimut flessibile (config_hub_telescopeonring_alhub_fit_ground_aztel.yaml)

Banda del fit 1-150 Hz, errore logaritmico complesso, perdita soft_l1, f_scale 0.1; anello di corrente 1500 Hz e ritardo equivalente; inerzie fisse al CAD. Modello: `src/gimbal_mdof/modi_forcella.py`.

## Polare: masse sullo statore

Misure: `BodeMeasurement_Polar_14vit_600_20261028_1perc.csv`, `BodeMeasurement_Polar_14vit_30_150_20261028.csv`, `BodeMeasurement_Polar_14vit_78_20261028_2perc_12mm_tightenedlaser_newsupport_4A_80degs_16screws.csv`

| massa | f [Hz] | zeta | a = kappa/J [1/(kg m^2)] | J [kg m^2] | kappa | k [N m/rad] |
|---|---|---|---|---|---|---|
| 1 | 8.84 | 0.1680 | +0.315 | 3.17 | +1 | 9.77e+03 |
| 2 | 53.37 | 0.0369 | +0.220 | 4.54 | +1 | 5.1e+05 |
| 3 | 57.47 | 0.0411 | -0.872 | 1.15 | -1 | 1.49e+05 |
| 4 | 100.24 | 0.0348 | +0.702 | 1.42 | +1 | 5.65e+05 |

- inerzia rigida dell'asse (CAD): 0.8998 kg m^2 (1/J = 1.111)
- telescopio - mozzo: k_b = 7.964e+05 N m/rad (H7 1.16e+05), c_b = 717.1 N m s/rad
- cuscinetti sulla struttura dello statore: k_g = 136.9 N m/rad, eta = 2.7
- ritardo equivalente: 0.683 ms
- costo 12.507; con i cuscinetti a terra (come H7) 13.386
- parametri al limite: nessuno
- poli del modello (Hz, zeta): 2.0 (0.001), 8.9 (0.167), 53.4 (0.037), 57.5 (0.041), 100.2 (0.035)
- zeri del modello (Hz, zeta; zeta < 0 = fase non minima): 7.8 (+0.148), 52.3 (+0.036), 74.2 (-0.224), 74.5 (+0.327), 182.9 (+0.438)

## Azimut: load flessibile

Misura: `BodeMeasurement_Azimuth_14vit_600_20261028.csv`

- base della forcella J_L = 0.4451, sommita' J_A = 0.8604 kg m^2 (alpha = 0.659); telescopio J_P = 0.8578 kg m^2 (CAD)
- base - sommita': k_LA = 4.179e+04 N m/rad, c_LA = 0.001 N m s/rad
- sommita' - telescopio: k_3 = 2.861e+04 N m/rad (H7 1.16e+05), c_3 = 278.3 N m s/rad
- albero motore - base: k_2 = 1.884e+04 N m/rad (H7 1.309e+04), c_2 = 0.001
- cuscinetti: k_g = 914.4 N m/rad, eta = 0.645; ritardo 0.453 ms
- costo 6.112; costi delle partenze: 6.112, 6.112, 6.112, 6.703, 6.703, 6.703, 6.703, 6.768, 6.770
- parametri al limite: `c_LA`, `c_2`
- modi (Hz, zeta, forma e/m/L/A/P): 3.2 (0.000: +0.97 +0.97 +0.97 +0.99 +1.00); 55.8 (0.052: +1.00 +1.00 +1.00 -0.28 -0.23)

## Qualita' del fit (NRMSE %, nella banda del fit; RMS degli errori)

| modello | misura | tutta | 1-5 Hz | 5-100 Hz | >100 Hz | RMS dB | RMS deg |
|---|---|---|---|---|---|---|---|
| H7 + 4 stator masses | elevation, 1 % grid (1-600 Hz) | 69.0 | 61.0 | 66.0 | 71.5 | 1.52 | 11.5 |
| H7 + 4 stator masses | elevation, 1 % grid (30-150 Hz) | 65.9 | nan | 64.2 | 81.8 | 1.87 | 6.9 |
| H7 + 4 stator masses | elevation, 16 screws, 4 A (30-300 Hz) | 49.4 | nan | 46.8 | 73.6 | 1.29 | 11.6 |
| H7 + stator masses, bearings to ground | elevation, 1 % grid (1-600 Hz) | 66.2 | 57.8 | 62.2 | 71.9 | 1.62 | 12.5 |
| H7 + stator masses, bearings to ground | elevation, 1 % grid (30-150 Hz) | 66.0 | nan | 64.4 | 82.1 | 1.87 | 6.8 |
| H7 + stator masses, bearings to ground | elevation, 16 screws, 4 A (30-300 Hz) | 49.8 | nan | 47.2 | 73.0 | 1.3 | 11.5 |
| H7 configuration | elevation, 1 % grid (1-600 Hz) | 37.9 | 30.1 | 24.8 | -181.6 | 9.06 | 52.7 |
| H7 configuration | elevation, 1 % grid (30-150 Hz) | -13.0 | nan | 10.6 | -167.1 | 15.4 | 86.9 |
| H7 configuration | elevation, 16 screws, 4 A (30-300 Hz) | -11.7 | nan | 8.3 | -173.4 | 15.62 | 88.8 |
| flexible load | azimuth, 3 % grid (1-600 Hz) | 69.6 | 68.2 | 62.0 | 6.7 | 2.14 | 8.6 |
| H7 configuration | azimuth, 3 % grid (1-600 Hz) | -21.1 | -44.9 | -16.0 | -556.8 | 4.41 | 30.6 |

Figure: `figures/bode_polar_fit.png`, `figures/bode_azimuth_fit.png`, `figures/polar_stator_masses.png`, `figures/fit_errors.png`
Tabelle: `tables/polar_parameters.csv`, `tables/azimuth_parameters.csv`, `tables/fit_quality.csv`, `tables/features.csv`, `tables/starts.csv`
