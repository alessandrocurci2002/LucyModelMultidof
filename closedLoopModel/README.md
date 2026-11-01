# Modello in anello chiuso di LUCY

Modello e studi di controllo del pointer: impianto meccanico + motori Kollmorgen TBM(S)-7646-A + drive Kollmorgen
AKD2G. La cartella contiene la libreria `lucy_cl` e gli script che la usano per:

- analizzare e validare l'anello chiuso sulle misure;
- identificare guadagni efficaci, plant e attrito;
- simulare l'inseguimento con osservatore di disturbo, EtherCAT e giroscopio MEMS.

Il modello meccanico viene da [`gimbal_mdof`](../README.md). Il diario dettagliato di decisioni e risultati è
[`STATO_LAVORI.md`](STATO_LAVORI.md).

## Stato

Aggiornato al 2026-11-01.

| attività | stato | risultati | report |
|---|---|---|---|
| Libreria e validazione di settembre | fatto | `results_validation_H6_zeta10/` | `report/` |
| Studio delle prove di ottobre, iterazioni I1–I5 | fatto | `results_VClosedLoop_H6_zeta10/`, `results_vcl_study/` | `report/`, `report_vcl/` |
| Plant pilotato dalla corrente, I6–I9b | fatto | `results_current_input*/`, `results_identifiability/` | `report_vcl/` |
| Identificazione sui doppietti V2 | fatto | `results_dataset_v2/`, `results_v2_identification/` | `report_v2/` |
| Osservatore di disturbo (Hilkert), continuo | fatto, solo simulazione | `results_dob_tracking/` | `report_dob/` |
| Stesso con EtherCAT 1 ms e ADIS16507 | fatto, solo simulazione | `results_dob_ethercat/` | `report_dob/` |
| Stimatore grey-box in stile `greyest` | fatto | `results_vtesi_greybox/`, `results_vcl_greybox*/` | nessuno; vedi `STATO_LAVORI.md` |
| **LQG alla Gawronski** | **audit fatto, sintesi e simulazione da fare** | nessuno | `LQG_GAWRONSKI_AUDIT.md` |

Nessun controllore di questa cartella risulta provato sull'hardware: osservatore, EtherCAT e giroscopio sono
simulazioni.

Nel codice non c'è alcuna sintesi LQG o LQR. Il prossimo passo, definito in `prompt.txt`, è sintetizzare e simulare il
controllore, con accuratezza e jitter sulla stessa circonferenza dell'osservatore. La mappatura verso il drive è già
ricavata nell'audit:

- K_pv = k_D;
- K_pp e ω_i sono le radici di z² − (k_P/k_D)·z + k_I/k_D, con K_vff = 1;
- K_f per lo smorzamento attivo è calcolato nel master EtherCAT.

## Avvio rapido

Tutti i comandi partono da questa cartella. L'ambiente è quello del progetto padre: l'installazione è descritta nel
[README principale](../README.md#installazione).

```bash
cd sim/gimbal_multidof/closedLoopModel
../.venv/bin/python -m pytest tests -p no:cacheprovider                 # test
../.venv/bin/python -I main.py --config config_zeta10.yaml --skip-time  # solo analisi in frequenza
../.venv/bin/python -I dob_tracking_study.py --replot                   # metriche e figure dalle run salvate
```

Le misure stanno in `../../LucyIdentification`, fuori dal repository e non sotto git. Le run salvate in `runs/*.npz`
evitano di rifare le simulazioni più lunghe.

## Struttura

```
closedLoopModel/
├── lucy_cl/                 libreria
├── data/motor_TBM7646A.yaml datasheet del motore
├── *.py                     script di studio (tabella sotto)
├── config*.yaml             una configurazione per variante o iterazione
├── results_*/               uscite, una cartella per configurazione
├── report*/                 report LaTeX generati
├── tests/                   test
├── tools/diag_sticking.py   diagnosi di un problema risolto, eliminabile
├── STATO_LAVORI.md          diario del lavoro
├── LQG_GAWRONSKI_AUDIT.md   materiale per la sintesi LQG
└── prompt.txt               prossimi passi concordati
```

## Libreria `lucy_cl`

| modulo | contenuto |
|---|---|
| `plant.py` | `GimbalPlant`: modello concentrato non lineare da qualsiasi configurazione di `gimbal_mdof`, attrito di Dahl |
| `motor.py` | `MotorParams`, `Motor`: corrente q per fase rms, K_t·i − b·ω |
| `driver.py` | `DriveParams`, `Drive`: posizione → velocità → corrente, biquad, osservatore |
| `controller_tf.py` | funzioni di trasferimento razionali del controllore, verificate sulle tracce WorkBench |
| `closed_loop.py` | `ClosedLoop`: simulazione Radau, linearizzazione, FRF, margini |
| `data_files.py` | lettori di `.akd` ed export WorkBench |
| `validation.py` | sequenze misurate, stato iniziale dalla misura, simulazione di finestre |
| `identification.py` | stima dei guadagni efficaci e del livello d'attrito |
| `current_driven.py` | `CurrentDrivenPlant`: plant pilotato dalla corrente misurata, surrogato rigido |
| `plant_fit.py` | fit dei parametri del plant su una FRF misurata |
| `greyest.py` | stima grey-box di modelli ODE lineari, in stile `greyest` di MATLAB |
| `gimbal_grey.py` | adattatore plant → `greyest` |
| `disturbance_observer.py` | `GyroObserverLoop`: osservatore a stato ridotto di Hilkert, notch opzionale |
| `ethercat.py` | `SampledObserverLoop`: anello a dati campionati con master a 1 ms |
| `gyro.py` | `MemsGyro`: ADIS16507 con rumore da varianza di Allan, `allan_deviation` |
| `output.py` | figure e tabelle di `main.py` |

Il segnale per asse: il driver legge il comando di posizione e la posizione e velocità dell'encoder lato carico, e dà
la tensione di fase. Il motore ne ricava la corrente e la coppia sul nodo motore. Il plant restituisce le posizioni.
Non ci sono loop algebrici.

Per cambiare il modello meccanico basta cambiare `plant.config` nel file di configurazione, usando un qualsiasi config
di `gimbal_mdof`.

## Script, configurazioni e risultati

| script | scopo | configurazione → risultati |
|---|---|---|
| `main.py` | componenti, driver contro WorkBench, FRF, margini, scenari | `config.yaml` → `results_H6/` (**mai eseguito**) |
| | | `config_zeta10.yaml` → `results_H6_zeta10/` |
| | | `config_vcl.yaml` → `results_VClosedLoop_H6_zeta10/` |
| `validate.py` | stima dei guadagni su una sequenza, verifica sulle altre | `config_validation.yaml` → `results_validation_H6_zeta10/` |
| `compare.py` | step di ottobre con i parametri dei driver, varianti di lettura | `config_compare_vcl.yaml` → `results_VClosedLoop_H6_zeta10/time_comparison/` |
| `vcl_study.py` | studio di ottobre, iterazioni I1–I5; scrive `config_vcl_fit.yaml` | `config_vcl_study.yaml` → `results_vcl_study/` |
| `current_input_study.py` | plant pilotato dalla corrente, screening e fit | `config_current_input.yaml`, `_I6`, `_I7`, `_I8`, `_I9a`, `_I9b` → `results_current_input*/` |
| `identifiability.py` | Jacobiano, condizionamento, errori standard | run fissate in `current_input` → `results_identifiability/` |
| `dataset_assessment.py` | contenuto informativo di un set di sequenze | I9a → `results_dataset_v2/` |
| `v2_identification.py` | scelta della struttura con cross-validation, fit e validazione; riprende se interrotto | `config_v2_identification.yaml` → `results_v2_identification/` |
| `cable_torque_analysis.py` | coppia resistente contro angolo e storia | stessa configurazione → `results_v2_identification/cables/` |
| `dob_tracking_study.py` | inseguimento, encoder contro osservatore; `--replot` | `config_dob_tracking.yaml` → `results_dob_tracking/` |
| `dob_ethercat_study.py` | stesso con EtherCAT e rumore del giroscopio; `--replot` | `config_dob_ethercat.yaml` → `results_dob_ethercat/` |
| `greybox_fit.py` | stima grey-box di H7 sulle FRF; scrive `*_fit.yaml` | `config_vtesi_greybox.yaml` → `results_vtesi_greybox/` |
| | | `config_vcl_greybox.yaml` → `results_vcl_greybox/` |
| | | `config_vcl_greybox_complex.yaml` → `results_vcl_greybox_complex/` |

Altri esempi:

```bash
../.venv/bin/python -I validate.py --config config_validation.yaml
../.venv/bin/python -I validate.py --config config_validation.yaml --skip-identification
../.venv/bin/python -I current_input_study.py --config config_current_input_I9b.yaml
../.venv/bin/python -I v2_identification.py --config config_v2_identification.yaml   # circa 150 min, riprende
../.venv/bin/python -I dob_ethercat_study.py --config config_dob_ethercat.yaml       # circa 70 min
../.venv/bin/python -I greybox_fit.py --config config_vcl_greybox_complex.yaml --axes azimuth
```

### Dipendenze fra i risultati

- `config_vcl_fit.yaml` è scritto da `vcl_study.py` e usato da `dob_tracking_study.py`. La sua `results_dir`,
  `results_vcl_fit/`, non esiste: nessuna analisi `main.py` l'ha usata.
- `config_vcl_greybox_fit.yaml`, `config_vcl_greybox_complex_fit.yaml` e `config_vtesi_greybox_fit.yaml` sono scritti da
  `greybox_fit.py`; le loro `results_dir` `*_closed_loop/` non esistono.
- `dob_tracking_study.py` usa i guadagni di `results_vcl_study/` e l'attrito di `results_v2_identification/`.
- `results_v2_identification/summary.json` è anche ingresso di `analisi_forcella_lugre.py` e degli altri script Dahl
  nella cartella padre.
- `lucy_cl` è importato anche da `../analisi_anello_chiuso_ibrido.py`, che lo aggiunge a `sys.path`: non è un pacchetto
  installato.

## Report

Testo in inglese. Gli script leggono i `summary.json`, copiano o generano le figure e compilano con `latexmk`. Tutti
accettano `--no-latex`.

| report | pagine | contenuto | rigenerazione |
|---|---|---|---|
| `report/closed_loop_report.pdf` | 13 | modello contro misure, settembre e ottobre | `cd report && ../../.venv/bin/python -I build_closed_loop_report.py` |
| `report_vcl/vcl_report.pdf` | 13 | iterazioni I1–I9b, identificabilità | `cd report_vcl && ../../.venv/bin/python -I build_vcl_report.py` |
| `report_v2/v2_report.pdf` | 21 | identificazione V2, modello d'attrito, coppia resistente e cavi | `cd report_v2 && ../../.venv/bin/python -I build_v2_report.py` |
| `report_dob/dob_report.pdf` | 18 | osservatore, EtherCAT, giroscopio, reiezione dei disturbi, coppia richiesta | `cd report_dob && ../../.venv/bin/python -I build_dob_report.py` |

Le figure non sono sotto git. Dopo un clone vanno rigenerate con gli script di studio prima di ricompilare i report.

## Risultati principali

I dettagli sono nei report e in `STATO_LAVORI.md`.

**Guadagni dei drive.** I `.akd` non riproducono le sequenze. PL.KP va letto in 1/s, non come 2π·PL.KP, e l'anello di
posizione efficace è di circa 8 1/s su entrambi gli assi. In elevazione il file fornito non descrive la configurazione
attiva durante le prove.

**Identificazione su V2** (modello H6 completo, fit di posizione, mediana):

| caso | azimut | elevazione |
|---|---|---|
| identificazione V2 | 89.6 % | 87.8 % |
| validazione sui 4 CSV, θ₀ da V2 | 85.6 % | 49.7 % |
| validazione, θ₀ di sessione stimato su un altro file | 89.5 % | 85.4 % |

θ₀ è una proprietà della sessione (cavi o zero dell'encoder), da ricalibrare a ogni prova.

**Inseguimento** (simulazione, circonferenza di Attorrese: raggio 10 m, quota 5 m, centro a 100 m, 0.2 rad/s; un'orbita
con feedforward di velocità, EtherCAT 1 ms e ADIS16507):

| caso | E rms | E max | jitter > 1 Hz | jitter > 10 Hz |
|---|---|---|---|---|
| attuale, encoder | 26.6 mdeg | 131 mdeg | 99 µrad | 1.9 µrad |
| osservatore 20 Hz | 2.7 mdeg | 8.0 mdeg | 9.7 µrad | 2.6 µrad |

Il residuo con l'osservatore è la torsione della catena azimutale fra motore/encoder e forcella
(k_ψ2 = 13 091 N m/rad, cioè 4.4 mdeg per N m), che le maglie sugli encoder non vedono. È la motivazione per un
controllore che contenga il modello torsionale.

## Test

```bash
../.venv/bin/python -m pytest tests -p no:cacheprovider
```

25 test, tutti superati al 2026-11-01:

- `test_closed_loop.py` (11): componenti contro i riferimenti e coerenza dell'interconnessione;
- `test_greyest.py` (14): stimatore grey-box su sistemi con parametri noti.

`disturbance_observer.py`, `ethercat.py`, `gyro.py`, `current_driven.py` e `identification.py` non hanno test diretti.

## Ipotesi principali

- **Driver**
  - PI in serie, con le convenzioni verificate sulla traccia Controller di WorkBench.
  - PI di corrente con lo zero su R/L.
  - Observer modellato come passa-basso di banda OBS.BW.
  - Anti-windup a back-calculation: lo schema dell'AKD non è documentato.
  - Modello continuo, senza campionamento (a 1 ms solo in `ethercat.py`).
  - Feedforward di velocità `kvff` assente per default: le prove indicano che non è attivo. Gli studi di inseguimento
    lo attivano esplicitamente.
- **Motore**: commutazione ideale, modello per fase rms con K_t = √3·K_b, senza cogging. K_t = 0.258 N m/A_rms a
  freddo, 0.235 a caldo, ±10 %.
- **Plant**
  - Telescopio bilanciato (g = 0). Le coppie di mantenimento misurate non mostrano lo sbilanciamento di 11 N m previsto
    dal CAD.
  - Inerzia del rotore aggiunta al nodo motore.
  - Attrito ai cuscinetti nel tempo: elemento di Dahl, con |v| regolarizzato come √(v²+ε²) − ε, equivalente al vincolo
    a terra k(1+jη) di H6 all'ampiezza delle prove in frequenza. Da `current_driven.py` in poi anche curva di Stribeck.
- **H6 come identificato** (ζ₁ = 0.01, segnaposto Q8)
  - Il loop di elevazione risulta instabile a ~306 Hz, per il modo hub–telescopio su k_b; la misura invece è stabile.
  - Le simulazioni nel tempo usano quindi la variante con ζ₁ = 0.10.
- **Sequenze temporali**
  - La colonna AMPS vale 2 × la corrente in A rms: lo indicano saturazione e bilancio inerziale, non un dato del PLC.
  - La configurazione dei driver durante le prove non è disponibile, quindi i guadagni efficaci vengono stimati dai
    dati.
- **Osservatore e giroscopio**
  - Bias del giroscopio modellato come Gauss–Markov, con τ del minimo di Allan supposto 50 s.
  - Il modo strutturale del modello è a ~216 Hz, quello misurato a ~95 Hz: nel sistema reale il limite di banda sarà
    più basso di quello simulato.

## Limiti noti

- La corrente non è riprodotta bene quanto la posizione: oscillazione a ~8–10 Hz e ripple a ~50–60 Hz durante
  l'accelerazione.
- Il moto dell'azimut sposta l'elevazione ferma di 5–20 mdeg nella misura, circa 0.2 mdeg nel modello.
- Movimenti da 0.1–0.2° (pre-scorrimento) e oltre 13° restano fuori dal buon fit: un solo stato di Dahl non basta.
- Le stime grey-box di H7 in polare sono degeneri o non utilizzabili per k_b: la topologia non produce il polo a
  52.6 Hz sotto lo zero a 66 Hz. Il modello con le masse sullo statore, nella cartella padre, risolve questo punto.
- Gli studi in anello chiuso di questa cartella usano H6 e il Dahl di V2. Non usano ancora la struttura modale e la
  stima ibrida della cartella padre.
- `PlantGreyModel.damping_penalty` conta i modi sovrasmorzati come ζ = 1, senza gradiente: le partenze polari con
  anello = 1 si fermano subito.
