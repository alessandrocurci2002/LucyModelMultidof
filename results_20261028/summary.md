# Misure del 28/10/2026 e rigidezza anello-telescopio (config_hub_telescopeonring_alhub_fit_ground_aztel.yaml)

Polare con la lente a 12 mm: `BodeMeasurement_Polar_14vit_600_20261028_1perc.csv`, `BodeMeasurement_Polar_14vit_600_20261028.csv`, `BodeMeasurement_Polar_14vit_30_150_20261028.csv`; azimut: `BodeMeasurement_Azimuth_14vit_600_20261028.csv`. Riferimento: VTesi del 23/10.

- risonanza polare misurata (griglia all'1 %): **56.33 Hz** (-45.3 dB); sulle altre griglie 57.44, 56.75 Hz
- antirisonanza subito sotto: 49.4 Hz
- I_yy telescopio (CAD, fissa) = 0.8662 kg m^2
- configurazione: k_b = 1.16e+05 N m/rad -> zero 58.2 Hz, polo 294.4 Hz
- **k_b per il polo a 56.33 Hz: 4052 N m/rad** -> zero 10.89 Hz
- inerzia dal lato motore contro il telescopio: J = I_yy / ((f_p/f_z)^2 - 1) = 0.0336 kg m^2
- k_b che metterebbe invece lo zero a 49.4 Hz: I_yy (2 pi f)^2 = 8.355e+04 N m/rad

Fit NRMSE (%), bande `all`, `5-100 Hz`:

- elevation, lens 12 mm, 1 % grid, 28 Oct: k_b = 1.16e+05 -463.9 / 24.8; k_b = 4.05e+03 -529.1 / -1276.8
- elevation, lens 12 mm, 3 % grid, 28 Oct: k_b = 1.16e+05 -448.9 / 24.9; k_b = 4.05e+03 -525.5 / -1258.0
- elevation, lens 12 mm, 30-150 Hz, 28 Oct: k_b = 1.16e+05 -13.0 / 10.6; k_b = 4.05e+03 -2792.7 / -2945.5
- elevation, 14 screws, 23 Oct: k_b = 1.16e+05 -488.5 / 19.7; k_b = 4.05e+03 -609.1 / -1289.1
- azimuth, 3 % grid, 28 Oct: configuration -23.9 / -16.0
- azimuth, 14 screws, 23 Oct: configuration 17.8 / -13.4

Figure: `figures/bode_polar_20261028.png`, `figures/bode_azimuth_20261028.png`, `figures/bode_measured_20261028.png`, `figures/k_b_sweep.png`
