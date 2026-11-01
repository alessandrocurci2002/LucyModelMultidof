# Telescopio appeso al load di azimut (config_hub_telescopeonring_alhub_fit_ground_aztel.yaml)

Modello concentrato linearizzato a theta0 = 0, canale coppia -> encoder di azimut. Base: config_hub_telescopeonring_alhub_fit_ground.yaml.

- k (L - telescopio) = 1.16e+05 N m/rad (first try: k_b of the polar axis, 12 ring-telescope screws, 0.2 mm uniform gap (58.2 Hz, 19 Oct 2026))
- J_L (forcella + alberi polari + mozzo + anello) = 1.3055 kg m^2, J_telescopio(theta0 = 0) = 0.8578 kg m^2, k_psi,2 = 13091 N m/rad
- **coppia del telescopio: polo 75.34 Hz, zero 76.02 Hz** (zero sopra il polo del 0.9 %)
- forma del modo (azimut, telescopio = 1): psi_e -0.664, psi_m -0.661, psi_L -0.657, psi_p +1.000
- beta di azimut (zeta sul primo modo flessibile, ora il modo del telescopio): H6 4.11e-06, H7 4.23e-05 s: gli smorzatori di azimut crescono di 10.3 volte

| modello | poli [Hz] | zeri [Hz] |
|---|---|---|
| H6 | 774.05, 1846.74 | 12.38 |
| H7 | 75.34, 774.09, 1846.74 | 12.27, 76.02 |

Misure (20-150 Hz, prominenza 4 dB):

- azimuth, 14 screws, 23 Oct 2026, 3.5 A: risonanza 48.9 Hz, risonanza 55.7 Hz, antirisonanza 62.1 Hz, antirisonanza 64.0 Hz
- azimuth, 20 Oct 2026, 3.5 A: risonanza 44.1 Hz, risonanza 49.6 Hz, antirisonanza 59.1 Hz
- azimuth, polar at 0 deg, 3.5 A: risonanza 37.4 Hz, antirisonanza 45.3 Hz

Scansione di k (2e+04-4e+05 N m/rad, `tables/k_sweep.csv`): lo scarto fra zero e polo resta fra 0.3 e 5.8 % (0.9 % a k = 1.16e+05).

Figure: `figures/bode_azimuth_H6_H7.pdf`, `figures/k_sweep.pdf`.
