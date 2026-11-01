# Modi della forcella: stima modale ai minimi quadrati, configurazione del 2026-10-31

Copia di sicurezza del modello riportato in `fork_modes_report.pdf` (quinto passo). Le cartelle di lavoro
(`results_fork_modes_modal`, `report_fork_modes`) vengono riscritte a ogni nuova esecuzione: questa no.

## Contenuto

| file | cosa |
|---|---|
| `config_hub_telescopeonring_alhub_fit_ground_aztel.yaml` | configurazione H7 usata (CAD, materiali, telescopio, K_t = 0.258 N m/A_rms) |
| `fork_modes_report.pdf` | report (inglese, 15 pagine) |
| `results/` | `summary.json` (parametri, modi, qualita' del fit), tabelle CSV e figure |
| `code/forcella_modale.py` | modello modale e stima (copia di `src/gimbal_mdof/forcella_modale.py`) |
| `code/analisi_forcella_modale.py` | script che produce `results/` |

## Impostazioni della stima

- banda del fit 1-150 Hz; polare: tre misure nominali; azimut: `BodeMeasurement_Azimuth_14vit_600_20261028.csv`;
- formulazione modale con cuscinetti lineari isteretici k_g (1 + j eta_g); telescopio rigido sull'anello;
- minimi quadrati sulle parti reale e immaginaria (nessun logaritmo, perdita quadratica);
- smorzamento modale zeta_k <= 0.05;
- inerzie e rigidezze del CAD entro +-20 %;
- K_t entro la tolleranza di catalogo +-10 % (fattore g_Kt per asse);
- prima partenza: stima con CAD +-20 % e K_t nominale (`results_fork_modes_modal_cad20`).

## Risultati principali (banda 1-150 Hz)

| | fit (NRMSE) | errore rms |
|---|---|---|
| polare, 1-600 Hz | 59.6 % | 2.70 dB |
| polare, 30-150 Hz | 68.6 % | 2.42 dB |
| polare, 16 viti | 65.2 % | 2.09 dB |
| azimut | 74.9 % | 2.51 dB |

- polare: masse sullo statore a 56.7 e 99.0 Hz (zeta 0.033), cuscinetti k_g = 129 N m/rad, eta_g = 2.72;
  g_Kt = 1.007 (non separabile dall'inerzia mozzo + telescopio);
- azimut: fork mode a 56.2 Hz (zeta 0.028); g_Kt = 0.90 al limite inferiore; inerzia del carico e k_2 al limite +20 %.

## Riprodurre

Da `sim/gimbal_multidof`, con il codice del progetto allo stato del 2026-10-31:

    .venv/bin/python -I analisi_forcella_modale.py

richiede `results_fork_modes_modal_cad20/summary.json` come prima partenza; senza, la stima parte dalle partenze di
controllo e arriva allo stesso minimo.
