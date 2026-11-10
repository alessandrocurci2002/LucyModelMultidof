# Stato dei lavori: closedLoopModel (aggiornato il 2026-10-19)

## Fatto

**Libreria `lucy_cl`:**
- `plant`, `motor`, `driver`, `closed_loop`, `output`;
- `validation` (sequenze temporali) e `identification` (stima dei guadagni efficaci).

**Script:**
- `main.py`: analisi del closed loop;
- `validate.py`: validazione nel tempo.

I test sono in `tests/test_closed_loop.py` (11, tutti passano) e c'è il `README.md`.

**Solutore:**
- Radau con rtol 1e-4 e atol × 100;
- Jacobiano fornito, con passi fissi per gruppo di stati;
- anti-windup a back-calculation (T_AW = 1 ms; 0.2 ms per il loop di corrente);
- |v| della Dahl regolarizzato con √(v²+ε²) − ε;
- ogni finestra di simulazione gira in 2–4 s.

## Validazione sulle sequenze (`validate.py --config config_validation.yaml`)

**Sequenze usate** (decisione dell'utente, 2026-10-19):
- identificazione su `MPCLucy_pos_dist_dataset`;
- validazione su `MPCLucy_pos_dist_test` e `MPCLucy_scotched`;
- `Position_dist_FF*.csv` non usati.

**Unità della colonna AMPS:** 2 × I_rms. Lo indicano la saturazione a 36 (= 2 × IL.LIMITP) e il bilancio inerziale; l'utente ha accettato la stima.

**Guadagni dei driver:** stimati dai dati, per scelta dell'utente, perché i .akd non riproducono le sequenze. L'errore di inseguimento è ~2° a 15 °/s su entrambi gli assi: con i .akd sarebbe 0.24° in azimut e 6° in elevazione.

**Parametri stimati** per asse: Kpp, Kpv, Kiv e Tc (Coulomb della Dahl).

Esecuzioni:

| file | stato | note |
|---|---|---|
| `summary_run1` | superato | stato iniziale difettoso |
| `summary_run2` | superato | senza deformata statica |
| run 3 (`summary.md`) | **finale**, convergente | |

**Stima finale (run 3):**

| asse | Kpp [1/s] | Kpv [A/(rad/s)] | Kiv [Hz] | Tc [N m] |
|---|---|---|---|---|
| azimut | 7.43 | 198 | 0.148 | 0.30 |
| elevazione | 8.58 | 61.4 | 0.73 | 0.29 |

Fit NRMSE su dati non visti, mediana:

| | azimut | elevazione |
|---|---|---|
| errore di posizione, parametri stimati | 92–93 % | 82–87 % |
| errore di posizione, .akd | da −17 % a +3 % | da −26 % a −8 % |
| corrente, parametri stimati | 37–42 % | 25–31 % |

Asse fermo durante il moto dell'altro: misura 5–29 mdeg, modello 0.1–0.2 mdeg.

**Discrepanze trovate:**
- Kpp(azimut) ≈ 7.4 1/s ≈ 0.74 × PL.KP: il guadagno di posizione sembra in 1/s, non 2π·PL.KP.
- L'elevazione nelle sequenze aveva una configurazione diversa dal file .akd dell'11/09.
- Attrito:
  - Tc ≈ 0.3 N m, contro 1.6 / 0.66 N m della Dahl tarata sulle FRF;
  - la coppia di mantenimento a riposo supera Tc (stiction), cosa che la Dahl non rappresenta.
- Telescopio di fatto bilanciato: g = 0 è confermato dai dati.
- Accoppiamento: il moto dell'azimut sposta l'elevazione ferma di 5–20 mdeg, contro ~0.2 mdeg nel modello.
- Corrente: il modello non riproduce l'oscillazione a ~8–10 Hz né il ripple a ~50–60 Hz durante l'accelerazione.
- Kpv dell'elevazione è limitato dal modo spurio a 306 Hz di H6.

## Prove di ottobre (`LucyIdentification/VClosedLoop`, 2026-10-20)

**Configurazioni:**
- `config_vcl.yaml` → `results_VClosedLoop_H6_zeta10/`, con i Bode e i driver nuovi;
- `config_compare_vcl.yaml` + `compare.py` → `results_VClosedLoop_H6_zeta10/time_comparison/`, con gli step da 1° e 2°.

**Parametri dei file:**
- azimut uguale a settembre;
- elevazione VL 25/15 con notch a 50/96 Hz; PL.KP 0.4, KVFF 1.

**Bode:**
- l'elevazione è coerente con la misura precedente (12 viti);
- l'azimut coincide sopra 10 Hz, ma sotto 10 Hz il picco del vincolo a terra scende da ~4 a ~2 Hz: il comportamento a bassa frequenza non è ripetibile.
- La Dahl prende l'ampiezza dalle misure su cui H6 è stato tarato (`amplitude_measurement`).

**Step, fit dell'errore di posizione (mediana):**

| variante | azimut | elevazione |
|---|---|---|
| file + 2π·PL.KP | 10 % | −9 % |
| **file + PL.KP in 1/s** | **76 %** | −275 % |
| guadagni stimati a settembre (H) | 78 % | **91 %** |

**Conclusioni:**
- I parametri forniti riproducono l'azimut se PL.KP è letto in 1/s.
- L'elevazione si comporta come i guadagni stimati a settembre (Kpp 8.6, VL.KP 61, VL.KI 0.73), non come il file.
- Il feedforward (KVFF = 1 nei file) peggiora sempre il fit: nelle prove non è attivo.

## Report

Il report (inglese, 13 pagine) è `report/closed_loop_report.pdf`. Si rigenera con:

```bash
cd report && ../../.venv/bin/python -I build_closed_loop_report.py
```

Lo script legge i `summary.json` di `validate.py`, `compare.py` e `main.py` (config_vcl) e rigenera figure, tabelle e numeri.

## Studio VClosedLoop (`vcl_study.py`, 2026-10-20)

**Configurazione:** `config_vcl_study.yaml` → `results_vcl_study/`, più `config_vcl_fit.yaml` con il plant ottimizzato.

**Passi 1–2, controllori:** le funzioni di trasferimento ricavate dai file .akd coincidono con le tracce WorkBench a meno di 0.001 dB.

**Passi 3–5, plant ottimizzato sui Bode di ottobre:**
- elevazione: k_b 1.16e5 → 1.97e5 N m/rad; anello κ_r 1000 → 1.9; molla a terra 164/3.11 → 188/2.76;
- azimut: molla a terra 1418/0.32 → 149/1.76;
- smorzatori per elemento con un prior debole.
- Le risonanze misurate a 95 Hz (elevazione) e 470 Hz (azimut) non si raggiungono con questi parametri: manca inerzia dal lato motore/hub.

**Passo 6, iterazioni I1–I5** (fit dell'errore / della corrente, mediane):

| iterazione | azimut | elevazione | note |
|---|---|---|---|
| I1: 2π·PL.KP | −6 / −146 | −10 / −8 | |
| I2: PL.KP in 1/s | 71 / −5 | −278 / −25 | |
| I3: guadagni identificati | 84 / 5 | 43 / −198 | |
| I4: + Tc | 83 / 45 | 42 / −15 | minimo locale in elevazione |
| I5: griglia con Tc misurato + stima congiunta | 83 / 45 | **81 / 68 (adeguato)** | |

**Valori finali (I5):**
- azimut: Kpp 6.98, VL.KP 379, VL.KI 0.26, Tc 0.17;
- elevazione: Kpp 6.78, VL.KP 54, VL.KI 2.9, Tc 0.39.

La corrente dell'azimut (45 %) resta sotto la soglia del 50 %. L'oscillazione misurata a 8–10 Hz non dipende da c_ψ2 né da VL.KP/VL.KI (verificato: 44.7–46.7 %), quindi è una caratteristica del plant che il modello non ha (per esempio la cedevolezza della base o dello statore). Va deciso con l'utente.

Lo sbilanciamento del telescopio non spiega la bassa frequenza del Bode di elevazione: servirebbe un offset di 1.25–1.95 m; il momento statico misurato è 0.25 N m. Dettagli in `results_vcl_study/imbalance/`.

## Studio con la corrente misurata come ingresso (`current_input_study.py`, 2026-10-20)

La corrente misurata pilota il plant (senza driver); l'unico indice è il fit NRMSE della posizione dell'encoder. Identificazione sugli step da 1°, validazione sugli step da 2°, first guess I5.

Aggiunte al modello:
- spunto dell'attrito: Dahl con curva di Stribeck LuGre (T_s, v_s);
- scala coppia/corrente per asse;
- ingresso filtrato a 150 Hz;
- stato iniziale con la gravità.

Risultati del fit di posizione in validazione (mediana; soglia 90 %):

| run | config | risultati | azimut | elevazione |
|---|---|---|---|---|
| primo | `config_current_input.yaml` | `results_current_input/` | **94.5 % (adeguato)** | 81 % |
| I6, sbilanciamento libero | `config_current_input_I6.yaml` | `results_current_input_I6/` | | 86.8 % |
| I7, + viscoso del cuscinetto | `config_current_input_I7.yaml` | `results_current_input_I7/` | | 86.2 % |

Dettagli per run:
- **Primo run.** Parametri azimut: scala 0.91, σ₀ 203, T_c 0.23, T_s/T_c 1.5. Lo screening mostra che k_b, l'anello e gli smorzatori non influenzano la posizione nel tempo.
- **I6.** d_x = −0.27 mm: lo sbilanciamento non spiega l'asimmetria.
- **I7.** b_θ,h = 0.04 N m s/rad. L'identificazione arriva al 96 %, ma i parametri stimati sugli step da 1° non generalizzano agli step da 2°: il modello arriva corto.

L'elevazione non è adeguata; i prossimi passi sono da decidere con l'utente.

**I8**, su richiesta dell'utente (`config_current_input_I8.yaml` → `results_current_input_I8/`):
- **Azimut.** Isteresi a due livelli (relè): H → 0.006 N m, nessun effetto; validazione 94.4 %.
- **Elevazione.** Coppia a·sin θ + b·cos θ + c al posto della molla torsionale (σ₀ reso rigido):
  - a = −18.4 N m, cioè una molla assoluta di 18 N m/rad; b + c = 0.21 N m;
  - identificazione 96 %, ma validazione **62.6 %**;
  - a 2° la coppia di posizione supera lo spunto e il modello scivola prima del movimento, mentre l'asse reale resta fermo.
- **Conclusione.** L'offset a regime non dipende dalla posizione assoluta né da un'isteresi a due livelli. Dipende dalla storia (pre-scorrimento che si azzera dopo lo strisciamento), come nel modello I7, che resta il migliore (86 %).

**I9**, molla dei cavi −K_c(θ − θ₀) con θ₀ incognito (ipotesi dell'utente), attrito di I7:
- **I9a, θ₀ unico per asse** (`results_current_input_I9a/`): validazione azimut 78.7 %, elevazione 71.5 %. Peggiora: la molla con punto neutro fisso non estrapola agli step da 2°.
- **I9b, θ₀ per sequenza** (`results_current_input_I9b/`), stimato da solo sulle sequenze di validazione dopo la correzione con griglia: validazione azimut **94.9 %**, elevazione **87.0 %**.
  - azimut: K_c 14.8 N m/rad, θ₀ tra 0.5 e 0.9° (coerente tra le prove);
  - elevazione: K_c 0.83 N m/rad, θ₀ tra 3.3 e 5.7°, cioè in pratica una coppia costante di 0.05–0.08 N m;
  - guadagno minimo rispetto a I7 e al primo run.

## Report dello studio VClosedLoop

Il report (inglese, 13 pagine) è `report_vcl/vcl_report.pdf`. Copre:
- parte A: closed loop, I1–I5;
- parte B: plant pilotato dalla corrente, dal primo run a I9b;
- ruolo dei parametri;
- analisi di identificabilità;
- prossimi passi.

Si rigenera con:

```bash
cd report_vcl && ../../.venv/bin/python -I build_vcl_report.py
```

**Prossimo passo concordato:** l'utente fornirà un nuovo dataset più vario. Poi si fa cross-validation per sequenza e si sceglie la struttura del modello in base al punteggio di cross-validation.

## Nuovo dataset V2 (`LucyIdentification/VClosedLoop/V2`, 2026-10-20)

16 file di doublet (+A, −2A, +A), un asse per file:
- ampiezze da 0.1 a 15° (azimut anche un ciclo fino a ±15°);
- 10 deg/s, pause di circa 35 s;
- finestre utilizzabili: azimut 26, elevazione 24.

Caratterizzazione in `results_dataset_v2/characterisation.json`.

Anomalie:
- `azimut_1deg_15deg_v2` contiene un solo movimento;
- `azimut_8deg_minus8deg` parte con un movimento in corso, già escluso dalle finestre.

**Contenuto informativo** (`dataset_assessment.py` → `results_dataset_v2/information.md`): modello I9a, errori standard allo stesso rumore, vecchio set da 1° contro V2.

| asse | condizionamento (vecchio → V2) | errori standard in V2 | correlazioni in V2 | fit di I9a su V2 |
|---|---|---|---|---|
| azimut | 17.5 → 66 | K_c 4 %, T_c 34 %, T_s/T_c 27 %, scala 13 %; σ₀ 127 % | T_c–T_s/T_c −0.95 | −38.7 % |
| elevazione | 215 → 14.5 | K_c 7 %, θ₀ 0.19°, T_c 34 %; v_s 150 %, b 86 % | nessuna sopra 0.9 | −2.5 % |

Sul vecchio set in elevazione le correlazioni erano T_c–b −0.99 e K_c–θ₀ −1.00.

I9a non estrapola ai movimenti grandi: la nuova identificazione deve ripartire da una ricerca a griglia.

**Proposta all'utente (in attesa di conferma):**
- identificazione su V2 (circa 49 finestre, senza `1deg_15deg_v2`), con cross-validation lasciando fuori un file alla volta per scegliere la struttura;
- validazione finale sui 4 CSV di VClosedLoop (26 finestre);
- fissare v_s e σ₀ (quest'ultimo eventualmente dal Bode), forse anche b;
- θ₀ per sequenza, stimato con una regola fissata prima di guardare la validazione.

## Identificazione su V2 (`v2_identification.py`, completata il 2026-10-21)

Lo schema di identificazione su V2 e validazione sui 4 CSV di VClosedLoop è stato confermato dall'utente.

**Configurazione:** `config_v2_identification.yaml` → `results_v2_identification/` (`run.log`, `summary.json`, `summary.md`, figure).

**Ripresa:** il run salva `summary.json` dopo ogni struttura. Rilanciando lo stesso comando riparte dalle strutture mancanti:

```bash
../.venv/bin/python -I v2_identification.py --config config_v2_identification.yaml
```

**Passi:**
0. Analisi senza modello delle pause e della crociera. Risultati preliminari:
   - l'attrito a 10 deg/s vale circa 0.37–0.38 N m su entrambi gli assi;
   - in azimut non c'è coppia di posizione;
   - in elevazione la coppia di posizione è circa 1.3–1.5 N m/rad;
   - la coppia di mantenimento dipende dalla direzione dell'ultimo movimento ed è ±0.2–0.35 N m.
   Il K_c di I9a (12–16 N m/rad) assorbiva in realtà la memoria dell'attrito.
1. Surrogato rigido (`CurrentDrivenPlant.simulate_window_rigid`, nuovo): ogni asse ruota come un corpo rigido. Coincide con il modello completo entro 0.1–0.7 punti di fit ed è 10–20 volte più veloce.
2. Strutture confrontate con cross-validation lasciando fuori un file alla volta (surrogato rigido):
   - molla sì/no;
   - poi σ₀, v_s e b fissati al nominale ×1/3, ×1, ×3;
   - poi le fissature accettate insieme.

   Regola di scelta: la mediana di CV più alta; a parità entro 1 punto, la struttura con meno parametri liberi. Sui fit si usa una perdita di Cauchy per finestra.
3. Fit della struttura scelta con il modello H6 completo su tutto V2.
4. Validazione sui 4 CSV di VClosedLoop con il modello completo:
   - θ₀ preso da V2;
   - in secondo luogo θ₀ stimato per sequenza.

**Risultati.** Durata 148 min, più 7 min per i controlli di validazione aggiunti. Il primo run è salvato in `summary_run1.json`.

Strutture scelte:
- **azimut:** Dahl + Stribeck + molla, σ₀ fissato a 4513 N m/rad (3 volte il valore dal Bode);
- **elevazione:** Dahl + Stribeck + molla, v_s fissato a 0.015 rad/s.

La molla migliora la cross-validation in entrambi gli assi:

| asse | CV senza molla | CV con molla |
|---|---|---|
| azimut | 81 % | 89 % |
| elevazione | 37 % | 84–87 % |

La mediana di CV resta entro 2.5 punti da quella di identificazione: niente overfitting.

Parametri finali (modello H6 completo):

| asse | scale | T_c [N m] | T_s/T_c | v_s [rad/s] | σ₀ [N m/rad] | K_c [N m/rad] | θ₀ [°] | b [N m s/rad] |
|---|---|---|---|---|---|---|---|---|
| azimut | 0.96 | 0.23 | 7.9 | 0.0058 | 4513 (fisso) | 0.93 | −2.68 | — |
| elevazione | 1.07 | 0.41 | 1.14 | 0.015 (fisso) | 1703 | 1.66 | −2.16 | ≈ 0 (sul limite) |

Fit di posizione (mediana):

| caso | azimut | elevazione |
|---|---|---|
| identificazione V2 | 89.6 % | 87.8 % |
| validazione sui 4 CSV, θ₀ da V2 | 85.6 % | 49.7 % |
| validazione, θ₀ di sessione stimato su un solo altro file | 89.5 % | 85.4 % |
| validazione, θ₀ per sequenza | 90.2 % | 86.4 % |

Osservazioni:
- θ₀ in elevazione vale +2.2…+2.9° nella vecchia sessione e −2.2° in V2: θ₀ è una proprietà della sessione (cavi o zero dell'encoder).
- In azimut θ₀ per sequenza è incoerente tra i file (−0.5…+2.2°) e la crociera non mostra coppia di posizione. Il K_c dell'azimut probabilmente compensa altro.
- Limiti residui:
  - movimenti da 0.1–0.2° (pre-scorrimento), con fit tra −3 e 69 %;
  - movimenti oltre 13°, con fit 73–85 %, per la deriva in crociera.

**Prossimi passi possibili (da decidere con l'utente):**
- attrito multi-stato (Maxwell-slip/GMS) per il pre-scorrimento;
- protocollo di calibrazione di θ₀ a ogni sessione;
- report LaTeX.

## Coppia resistente e report V2 (2026-10-21)

**Analisi della coppia resistente:** `cable_torque_analysis.py` → `results_v2_identification/cables/` (summary.md/json e figure). Usa i file V2 nell'ordine di registrazione, ricavato dalle date di modifica dei file. Risultati:

- **Crociera selezionata sulla velocità dell'encoder.** La selezione sulla velocità comandata includeva l'accelerazione, perché l'encoder insegue il comando con un ritardo di circa 1.5°. Con la selezione corretta:
  - coppia di posizione lineare in **entrambi** gli assi: azimut 1.31 N m/rad, elevazione 1.38 N m/rad. La conclusione precedente "azimut senza coppia di posizione" era sbagliata;
  - in elevazione la curva si ripete tra i file;
  - attrito a 10 deg/s: elevazione piatto, 0.36 N m; azimut da 0.4 a 0.2 N m lungo lo sweep.
- **Coppia di mantenimento.** Cresce con l'ampiezza dell'ultimo movimento:
  - azimut da 0.1–0.2 a 0.35–0.45 N m;
  - elevazione da 0.37–0.44 a 0.45–0.54 N m.

  È un'isteresi lenta, distribuita su alcuni gradi, che il singolo elemento Dahl non rappresenta.
- **Ritorno a 0 dopo il doublet.**
  - Elevazione: +0.08–0.18 N m subito dopo, poi in circa 1 min torna al valore precedente (0.36–0.41 N m), senza deriva lungo la sessione.
  - Azimut: dipende dall'ampiezza dell'ultimo doublet (0.18 o 0.36 N m), è reversibile e non si accumula.
- **Rilassamento a posizione costante.**
  - Elevazione: fino a 0.29 N m, τ circa 36 s.
  - Azimut: fino a 0.09 N m, τ circa 10 s.
  - L'origine (integratore del drive sotto un conteggio oppure meccanismo viscoelastico) non è decidibile dai dati registrati.
- **Tra le due sessioni,** in elevazione, la coppia di posizione si sposta di circa 0.15 N m. Lo conferma sia l'analisi senza modello delle pause sia θ₀ stimato per file. In azimut non c'è spostamento.

**Report** (inglese, 21 pagine): `report_v2/v2_report.pdf`. Si rigenera con:

```bash
cd report_v2 && ../../.venv/bin/python -I build_v2_report.py
```

Contenuto: metodologia passo passo, sezione sul modello d'attrito, sezione sulla coppia resistente e sui cavi.

Modifiche a `v2_identification.py`:
- validazione con θ₀ di sessione stimato su un solo file;
- jackknife;
- `cruise_torques(by_encoder=...)`.

## Da fare

1. (fatto) Risultati della run 3 riportati all'utente.
2. Decidere con l'utente la convenzione di PL.KP: le prove di ottobre confermano 1/s per l'azimut. Se si adotta:
   - aggiornare `DriveParams` (Kpp = PL.KP);
   - lasciare `controller_workbench`, che deve riprodurre la traccia WorkBench, su 2π·PL.KP;
   - rifare l'analisi in frequenza di `main.py`.
3. Eseguire `main.py` su `config.yaml` (H6 così com'è: elevazione instabile a 306 Hz) e su `config_zeta10.yaml`.
4. Elevazione: verificare sul drive, durante le prove, i valori effettivi di PL.KP, VL.KP e VL.KI. Il programma del PLC (MPC) potrebbe sovrascriverli a runtime.

`tools/diag_sticking.py` (diagnosi del problema di integrazione, ora risolto) si può eliminare.

## Osservatore di disturbo con giroscopio (Hilkert) sull'inseguimento della circonferenza (2026-10-23)

**Richiesta dell'utente:** confrontare due strutture su una traiettoria di inseguimento, con metriche di errore di puntamento e jitter:
- il controllo attuale, con maglie sugli encoder;
- lo stesso più un osservatore di disturbo a stato ridotto in feedforward (Hilkert & Pautler 2011), con un giroscopio virtuale (senza rumore) sul telescopio.

Riferimento: la circonferenza di Attorrese (raggio 10 m, quota 5 m, centro a 100 m, 0.2 rad/s).

**Codice:**
- `lucy_cl/disturbance_observer.py`: `GyroObserverLoop` (sottoclasse di `ClosedLoop`), `ObserverParams`, `rigid_inertia`. Notch opzionale sul percorso del feedforward.
- `ClosedLoop` ha ora `torque_scale`, `cable_spring` e `signal_names`, con default neutri; i test passano.
- `dob_tracking_study.py` con `config_dob_tracking.yaml` → `results_dob_tracking/` (`runs/*.npz`, `summary.md/json`, figure). L'opzione `--replot` rifà solo metriche e figure.

**Modello:**
- plant `config_vcl_fit.yaml`, guadagni dei drive di I5;
- attrito, molla dei cavi e scala di coppia da V2;
- osservatore con J′ = rᵀMr e K_t′ nominale.

**Stabilità:**
- In elevazione l'osservatore senza notch è stabile fino a 7 Hz di banda; da 8 Hz destabilizza il modo strutturale del modello a circa 216 Hz (telescopio contro motore).
- Con un notch a 210 Hz sul feedforward (come nella figura 3 di Hilkert) è stabile fino a 40 Hz, anche con il notch spostato del ±10 %.
- Il modo misurato è a circa 95 Hz: nel sistema reale il limite sarà più basso.

**Prima indicazione (finestra corta):**
- senza feedforward di velocità domina l'errore di inseguimento v/K_pp, circa 137 mdeg;
- con il feedforward di velocità (KVFF = 1) l'errore va da 7 a 1.1 mdeg aggiungendo l'osservatore.

**Risultati** (`results_dob_tracking/summary.md`), un'orbita dopo 5 s di assestamento, metriche sull'asse ottico del telescopio:

| caso | E rms [mdeg] | E max [mdeg] | jitter > 1 Hz [µrad] | jitter > 10 Hz [µrad] |
|---|---|---|---|---|
| attuale (encoder) | 124.1 | 183.7 | 99.4 | 1.89 |
| + osservatore 20 Hz | 117.7 | 182.6 | 4.4 | 0.010 |
| attuale + FF di velocità | 26.6 | 130.9 | 99.4 | 1.89 |
| + osservatore 20 Hz | 2.70 | 7.8 | 4.4 | 0.010 |

- **Senza FF di velocità** domina l'errore di inseguimento v/K_pp (circa 117 mdeg in azimut), che l'osservatore per costruzione non cambia (principio di dualità).
- **Con il FF** domina il disturbo, e l'osservatore riduce l'errore di 10 volte e il massimo di 17.
- **Il jitter scende di circa 22 volte** in entrambi i casi.
- **Residuo con l'osservatore:** è la torsione della catena azimutale tra motore/encoder e forcella (k_ψ2 = 13 091 N m/rad, cioè 4.4 mdeg per N m). Le maglie sugli encoder non la vedono.
- **Sensibilità:**
  - con il giroscopio a 100 Hz non cambia nulla;
  - con l'attrito d'azimut alternativo (T_s/T_c 2.66) l'errore va da 11.8 a 1.7 mdeg e il jitter da 3.6 a 0.4 µrad. I valori assoluti dipendono dal picco di Stribeck dell'azimut, il guadagno relativo resta.
- **Banda dell'osservatore:** oltre 5 Hz l'errore non cambia; il jitter > 1 Hz migliora fino a 20 Hz, poi i ritorni sono decrescenti.

## Vincoli sui segnali: EtherCAT a 1 ms e giroscopio ADIS16507 (2026-10-23)

**Codice:**
- `lucy_cl/ethercat.py`: `SampledObserverLoop`, simulazione ibrida.
  - Ciclo del master a 1 ms, con il solutore Radau mantenuto tra i cicli.
  - Target di posizione interpolati dal drive (lineare oppure ZOH).
  - Correnti lette dal TxPDO; osservatore discreto nel master; feedforward applicato al ciclo successivo, con il notch discreto.
- `lucy_cl/gyro.py`: `MemsGyro`, che modella l'ADIS16507-1 sugli assi x e y.
  - Dinamica del secondo ordine a 480 Hz, 2000 SPS, latenza 0.5 ms.
  - Rumore con ARW 0.13°/√h e instabilità del bias 2.3°/h (bias di Gauss–Markov, τ del minimo supposto 50 s).
  - Rumore in uscita 0.070°/s RMS, contro 0.068 della scheda tecnica.
  - Funzione `allan_deviation`.
- `dob_ethercat_study.py` con `config_dob_ethercat.yaml` → `results_dob_ethercat/`. Il run dura 70 min; i casi con osservatore e rumore costano circa 1 h di CPU ciascuno. L'opzione `--replot` rifà solo metriche e figure.

**Risultati** (un'orbita, con il FF di velocità):

| caso | E rms | E max | jitter > 1 Hz | jitter > 10 Hz | jitter > 100 Hz |
|---|---|---|---|---|---|
| attuale (con FF) | 26.6 mdeg | 131 mdeg | 99 µrad | 1.9 µrad | 0.0005 µrad |
| osservatore 20 Hz | 2.7 mdeg | 8.0 mdeg | 9.7 µrad | 2.6 µrad | 0.14 µrad |

- Il rumore del giroscopio raddoppia il jitter sopra 1 Hz rispetto al giroscopio ideale (da 5.0 a 9.7 µrad) e domina sopra 10 Hz.
- La banda migliore è 20 Hz: a 40 Hz il rumore cresce più del guadagno, a 5–10 Hz il disturbo resta.
- EtherCAT a 1 ms, con interpolazione lineare oppure ZOH, cambia poco l'errore di puntamento (2.70 → 2.72 mdeg, ZOH 2.90).
- La corrente di feedforward ha 2.9 A RMS in azimut, contro 2.55 A del giroscopio ideale.

**Report** (inglese, 14 pagine): `report_dob/dob_report.pdf`. Si rigenera con:

```bash
cd report_dob && ../../.venv/bin/python -I build_dob_report.py
```

Contenuto:
- schemi a blocchi in TikZ: struttura attuale, struttura proposta, osservatore di Hilkert, architettura EtherCAT/SPI;
- implementazione passo passo, in 7 passi;
- risultati nel continuo;
- focus su errore e jitter (attuale contro osservatore con giroscopio rumoroso): tabelle in mdeg/µrad e in gradi, jitter per banda, jitter cumulativo;
- limiti e prossimi passi.

Cita Hilkert & Pautler 2011, Hilkert 2008, Attorrese 2026, la scheda tecnica ADIS16507 e IEEE Std 952.

## Audit di Gawronski per l'LQG (2026-10-23)

`LQG_GAWRONSKI_AUDIT.md` raccoglie il materiale per sviluppare l'LQG di Gawronski (*Modeling and Control of Antennas and Telescopes*, 2008) mantenendo la struttura dell'AKD. Non è implementato.

Contenuto:
- mappa del libro con le pagine;
- equazioni (6.1–6.7, 8.1–8.27);
- mappatura esatta tra i guadagni LQR rigidi [k_I k_P k_D] e la cascata AKD: K_pv = k_D, K_pp e ω_i radici di z² − (k_P/k_D)·z + k_I/k_D (serve K_vff = 1);
- K_f calcolato nel master;
- numeri attuali di LUCY e procedura di sintesi in 8 passi;
- punti aperti.

Il prossimo passo (`prompt.txt`, seconda riga) è la sintesi e la simulazione, con le metriche di accuratezza e jitter sulla circonferenza di Attorrese.

**Aggiunte al report dell'osservatore** (2026-10-23, `report_dob`, ora 18 pagine):
- **§4.1, reiezione dei disturbi.** Definizione rigorosa di G_d = c(jωI − A)⁻¹b_d (coppia sul cuscinetto → asse ottico, attrito in strisciamento con σ₀ → 0), più:
  - modello rigido con i filtri del drive;
  - asintoti ω/(κK_t·k_I), ω²/(κK_t·k_I·ω_o), 1/(Jω²);
  - limite di cedevolezza 1/k_ψ2 = 76 µrad/(N m) in azimut;
  - tabella dei valori con il rapporto confrontato con il fattore di Hilkert;
  - verifica contro i massimi simulati (131 e 7.8 mdeg).
- **§5.3, coppia motore richiesta** (EtherCAT, FF). Attuale e osservatore con giroscopio ideale hanno la stessa coppia (0.65/0.63 N m RMS in azimut), ma l'osservatore la inverte prima. L'ADIS16507 aggiunge 0.33 N m RMS sopra 10 Hz e porta il picco a 2.69 N m. In elevazione il rumore fa da dither sul Dahl (attrito medio da 0.43 a 0.30 N m).

## Stima grey-box come `greyest` di MATLAB (2026-10-27)

`lucy_cl/greyest.py`, test in `tests/test_greyest.py` (12, passano). Non ancora usato sui dati di LUCY.

- **Modello** (`GreyModel`, come `idgrey`): `odefun(p, Ts) -> (A, B, C, D)`, continuo o discreto; parametri con valore,
  limiti, libero/fisso, ricerca sul logaritmo per i parametri positivi. In più `frfun` per modelli definiti solo in
  frequenza (vincolo a terra isteretico k(1 + j eta)).
- **Dati**: `FrequencyData` (come `idfrd`) o `TimeData` (come `iddata`, ZOH), anche più esperimenti.
- **Opzioni** (`GreyestOptions`): metodo (trf/lm/dogbox), filtro di pesatura (bande, vettore, inv, invsqrt), errore in
  frequenza complesso (come MATLAB) o logaritmico (dB + fase), peso sulle uscite, stato iniziale stimato, stabilità
  imposta, perdita robusta, regolarizzazione (Lambda, R, Nominal).
- **Report**: FitPercent (NRMSE), MSE, FPE, AIC/AICc/BIC, covarianza, std, correlazioni, condizionamento dello Jacobiano,
  parametri al limite; `compare`.
- **Covarianza**: predefinita la stima robusta (sandwich). Con rumore proporzionale alla risposta e errore complesso la
  formula classica di MATLAB sottostima la dispersione di circa 7 volte (verificato con Monte Carlo); con l'errore
  logaritmico le due coincidono.

## Stima grey-box di H7 sulle FRF VTesi (2026-10-27)

`greybox_fit.py` + `config_vtesi_greybox.yaml` -> `results_vtesi_greybox/` (summary.md/json, tabelle, Bode);
adattatore `lucy_cl/gimbal_grey.py` (verificato: riproduce `loop_frf` a 1e-16). Modello H7 a 9 GdL, errore logaritmico,
perdita soft_l1, 1-400 Hz, partenze multiple sulle molle libere. Vincoli dell'utente: inerzie e molle geometriche +-20 %,
k_theta,4, anello, k_psi,3 e vincoli a terra liberi, smorzamento libero. Il risultato e' scritto anche come
`config_vtesi_greybox_fit.yaml`.

- **Azimut** (fit 36 -> 51 %; 1-5 Hz -73 -> 75 %): k_psi,3 = 8.1e4 N m/rad (+-12 %), k_g = 1043, eta = 0.73.
  La coppia del telescopio resta stretta (polo ~57, zero ~59 Hz contro 55.7 / 62 misurati); 7 parametri su 13 al
  limite del +-20 %, J_e, J_m, k_psi,1, c_psi,3 non determinati.
- **Polare** (fit -511 -> 43 %): soluzione **degenere**. k_theta,4 = 16 N m/rad con c_theta,4 = 484 N m s/rad: il
  telescopio e' agganciato al mozzo da uno smorzatore invece che da una molla, cosi' sparisce l'antirisonanza (che la
  misura a 14 viti non ha) e resta una curva liscia che media le tre risonanze a 57 / 80 / 97 Hz. Tutti i parametri
  strutturali indeterminati (std > 50 %). Non fisico: conseguenza dello smorzamento libero con un modello che non
  contiene quelle risonanze.

## Stima grey-box di H7 sulle FRF del 20/10 con smorzamento fisico (2026-10-27)

`greybox_fit.py` (ex `vtesi_greybox_fit.py`, ora generico) + `config_vcl_greybox.yaml` -> `results_vcl_greybox/`;
configurazione in anello chiuso con la stima: `config_vcl_greybox_fit.yaml`. Misure `VClosedLoop/BodeMeasurement_*.csv`,
banda 1-150 Hz, vincoli come sopra piu' **zeta <= 0.1 su ogni modo dell'asse** (vincolo morbido, `penalty_fun` di
`greyest`, `PlantGreyModel.damping_penalty`); il fattore di perdita eta del vincolo a terra (attrito) resta libero.

- Primo tentativo con il vincolo solo sotto 300 Hz: l'ottimizzatore ha spostato modi appena sopra la soglia e li ha
  sovrasmorzati (300.0 / 304.5 / 343 Hz con zeta 1). Ora il vincolo vale su tutti i modi, separatamente per asse.
- **Azimut** (fit in banda -31 -> 54 %; 1-5 Hz -127 -> 85 %): k_psi,3 = 7.4e4 N m/rad (+-34 %), k_g = 215 N m/rad,
  eta = 1.22 (+-6 %); modi 55.0 Hz (zeta 0.004), 900 Hz (zeta 0.10, al limite), 1991 Hz. La coppia del telescopio
  resta stretta (polo/zero ~52/55 Hz contro 49.6/59 misurati); inerzie e molle degli alberi non determinate.
- **Polare** (fit in banda 46 -> 44 %): il telescopio diventa rigido sul mozzo (k_theta,4 = 1.1e7) e l'antirisonanza
  a 66 Hz sparisce; anche con l'errore complesso (prova non salvata) k_theta,4 -> 3.8e8 e fit 54 %. Con questa
  topologia il modello non produce il polo a 52.6 Hz sotto lo zero a 66 Hz della misura: per il fit e' meglio
  togliere l'antirisonanza che metterla nel posto sbagliato. Stima polare non utilizzabile per k_b.

## Stessa stima con l'errore complesso, come MATLAB greyest (2026-10-27)

`config_vcl_greybox_complex.yaml` -> `results_vcl_greybox_complex/` (e `config_vcl_greybox_complex_fit.yaml`):
errore complesso, nessuna pesatura, nessuna perdita robusta; resto invariato (1-150 Hz, +-20 %, zeta <= 0.1).

- **Azimut** (fit in banda -31 -> 58 %; 1-5 Hz 73 %, 5-100 Hz 38 %): k_psi,3 = 5.5e4 N m/rad (+-5 %), k_g = 147,
  eta = 1.83. Le inerzie vanno nell'angolo che allarga la coppia (J_L -20 %, J^ax -20 %, I_zz +20 %), ma la coppia
  resta a ~52 Hz e stretta, come previsto dalla catena a 4 GdL. L'ottimizzatore arriva a 3000 valutazioni
  scorrendo su una valle piatta (fit 57.4 -> 57.7 %): non convergenza formale.
- **Polare** (fit in banda 46 -> 54.5 %): telescopio rigido sul mozzo (k_theta,4 = 2e8), nessuna antirisonanza;
  come con l'errore logaritmico. Struttura insufficiente per il polo a 52.6 Hz sotto lo zero a 66 Hz.
- Partenze polari con anello = 1: si fermano subito (modi sovrasmorzati contati come zeta = 1, penalita' senza
  gradiente). Da migliorare in `PlantGreyModel.damping_penalty`.

## LQG di Gawronski (2026-11-01)

Sintesi in `lucy_cl/lqg.py`, studio in `../analisi_lqg_gawronski.py`, report in `../report_lqg_gawronski/`. Solo simulazione.

- Architettura PL del libro: plant = maglia di velocità chiusa, l'LQG sostituisce l'anello di posizione nel master a 1 kHz; stato $[\hat e_i, \hat e, \hat x_f]$. `dlqe` lavora su $[y, x_f]$: $y_i$ non è osservabile.
- Guadagni: azimut $k_p=14.8$, $k_i=24.1$, $q_f=3\cdot10^4$, $\rho_e=10^{-3}$; polare $k_p=129$, $k_i=123$, $q_f=0$.
- Punti aperti: rumore e quantizzazione dell'encoder, plant misurato con la configurazione reale del PLC, anti-windup, variante LP/LL, command shaping.

## Roadmap LQG-PI in stile DSN (definita dall'utente, 2026-11-01)

Obiettivo: usare il modello costruito come banco di prova per valutare l'efficacia di un controllore con LQG in posizione e PI in velocità in stile DSN (configurazione PL), confrontandolo con la cascata attuale. Tutte le fasi sono in simulazione.

1. **Taratura secondo Gawronski.** Per coerenza accademica, i pesi si tarano seguendo le indicazioni del libro (§8.2, §8.4). Riferimento: `LQG_GAWRONSKI_PL_TUNING_AUDIT.md`. Stato: da fare; la taratura attuale in `../results_lqg_gawronski/` non segue la procedura del libro su ordine dei passi e metriche.
2. **Inseguimento della circonferenza di Attorrese.** Stesso tipo di prova usata per validare il metodo di Hilkert per la reiezione dei disturbi (`dob_tracking_study.py`), con l'aggiunta di process noise e disturbi, per verificare se l'approccio LQG-PI è più robusto della cascata.
3. **Traiettorie reali del drone.** Confronto delle prestazioni di inseguimento fra schema di Gawronski e cascata attuale, su traiettorie con stimatori EKF2 e OpenVINS e groundtruth TagSLAM. Il riferimento al simulatore arriva dagli stimatori; l'errore di puntamento si calcola sulla groundtruth. Domanda: l'approccio LQG-PI riduce l'errore di puntamento?
   - Da reperire: i dati non sono ancora stati individuati in questa cartella.

### Decisioni per la taratura (2026-11-01)

Registro completo in `../report_lqg_tuning/`.

- Cuscinetti nel modello di `dlqr` e `dlqe` con la linearizzazione armonica. Decisione dell'utente: l'LQG reale può contenere solo un modello lineare, e la dissipazione non va sovrastimata ai modi.
- Elemento di Maxwell $K_b(s)=k_1s/(s+a)$, uno stato in più: non è altro che un sistema molla-smorzatore, molla $k_1$ in serie allo smorzatore $c=k_1/a$. Niente molla statica ($k_0=0$). Priorità alle risonanze, decisione dell'utente: peso 3 sulle bande di risonanza nel fit sul plant $y/u$ con la funzione descrittiva all'ampiezza della prova. Azimut $k_1=1268$ N m/rad, $c=244$ N m s/rad, $f_a=0.83$ Hz; polare $k_1=11824$ N m/rad, $c=34$ N m s/rad, $f_a=55.5$ Hz. Il viscoso $k_{eq}+c_{eq}$ è scartato: a 56 Hz in azimut porta lo smorzamento da 0.023 a 0.073.
- Il plant resta Dahl: Maxwell sta solo nel modello interno. Nel tempo perde la correlazione (anello aperto, fit mediano azimut da 42 a −86 %, polare da 13 a −52 %), perché durante il moto è viscoso e non satura a $T_c$. Nel passo 1 si tarano entrambi i modelli interni, Maxwell e senza cuscinetti, e si misurano i limiti elencati nel report.
- Punto di lavoro: ampiezza del cuscinetto nelle prove FRF. Da rivedere nella fase 2, confrontando con il progetto senza cuscinetti: sulla circonferenza senza rumore i cuscinetti strisciano per la maggior parte dell'orbita.
- Niente feedforward e niente CPP: coerenza con la formulazione del capitolo 8. Decisione dell'utente.
- Circonferenza: raggio 10 m, quota 5 m, centro a 100 m, 0.2 rad/s, come lo studio di Hilkert. Decisione dell'utente.

- **Passo 1, pesi flessibili (fatto, 2026-11-02):** §8.2 alla lettera, TF sul solo modello ridotto con Maxwell, $Q=V$, pesi bassi (decisioni dell'utente); $V$ sarà un parametro della taratura fine. Vettore: $q_i=0.1$, $q_p=0.3$, $q_f=0.1$ su ogni stato flessibile. I picchi (−59.5 dB az, −62.7 dB el) non sono troppo alti rispetto alla Fig. 8.5, quindi nessun peso alzato. Corretto un errore: $V$ entra in `dlqe` come $VT^2$ (`v_continua`), le conclusioni precedenti su instabilità e ciclo limite sono ritirate. Da tenere d'occhio: bias di 4 mdeg in elevazione su Dahl. Audit in `LQG_GAWRONSKI_PL_TUNING_AUDIT.md` §8. Prossimo: passo 2, $q_p$.

- **Passi 2-4, iterazione e $V$ (fatti, 2026-11-02):** azimut $q_i=41.9$, $q_p=19.4$, $q_f=0.1$ su tutti gli stati; elevazione $q_i=6.3\cdot10^4$, $q_p=2.5\cdot10^6$, $q_f=3\cdot10^4$ sul modo a 56 Hz e 0.1 sugli altri; $Q=V$. Sul modello l'obiettivo migliora molto in elevazione e poco in azimut. Sul plant Dahl l'azimut è peggiore della cascata (oscilla), l'elevazione molto migliore, ma con $k_p\approx1400$ e senza rumore. Prossimo: reiezione dei disturbi.

- **Confronto con le misure (fatto, 2026-11-02):** il banco riproduce la cascata misurata (errore rms 2 % sotto in azimut, 10 % in elevazione; la corrente meno). In azimut la perdita viene dallo stato di Maxwell nello stimatore: con modello interno senza cuscinetti l'LQG-PI batte la cascata ($k_p\approx25$, 255 mdeg contro 372). In elevazione il risultato è un limite superiore da banco senza rumore: 20 mdeg con $k_p=1425$, margini sul plant misurato bassi. Da decidere con l'utente: rivedere Maxwell nel modello interno. Script `analisi_lqg_confronto_misure.py`.

- **Regola (utente, 2026-11-02):** in simulazione i controllori si confrontano solo con il modello con Dahl, proxy della realtà, non con le prestazioni misurate.
- **Cuscinetti e libro:** i modi strutturali di Gawronski vengono dal FEM (§2.2, §8.1.4); l'attrito è trattato a parte (§11.2). Con modello interno senza cuscinetti e $q_i\approx0$, a pari errore l'LQG eguaglia la cascata e toglie la sovraelongazione (`analisi_lqg_equivalenza_cascata.py`). Da decidere: rifare i passi 2-4 con questo modello.

- **Taratura in vigore (2026-11-02, decisione dell'utente):** stato stimato con modello senza cuscinetti, simulatore Dahl. Azimut $q_i=0.25$, $q_p=645$, $q_f=0.1$; elevazione $q_i=1037$, $q_p=7.6\cdot10^5$, $q_f=10^4$ sul modo a 56 Hz; $Q=V$. Script con `--senza-cuscinetti`; risultati in `../results_lqg_tuning/no_bearings/`.
- **Fase 2, circonferenza (fatta):** `analisi_lqg_circonferenza.py`. Errore di puntamento 105 → 64 mdeg; con $k_i\approx0.5k_p$ 26 mdeg, con $q_p=10^5$ 17 mdeg. L'effetto dei disturbi si misura come differenza dal caso nominale, perché l'rms totale è dominato dall'inseguimento. In azimut il guadagno d'anello tarato è uguale alla cascata: $K_f$ e lo stimatore lento riassorbono $k_p$. La simulazione coincide con la $T_{wy}$ lineare. Accoppiamento fra gli assi 60-125 dB sotto il termine diretto. Da decidere: quali comandi fissano il limite di corrente (A1). Prossimo: fase 3, dati ancora da reperire.

- **Accoppiamento fra gli assi nel simulatore (fatto, 2026-11-02):** `simula_coppia` in `analisi_lqg_circonferenza.py` contiene coppie di Coriolis e giroscopiche, termini inerziali $H$ e $J_\psi(\theta)$ del telescopio, verificati sul modello simbolico di Lagrange. Sulla circonferenza del test effetto sotto 0.06 µrad. Su una circonferenza veloce e alta domina $J(\theta)$ (350 µrad rms). Il tuning è stabile fino a 80°, la variante $q_p=10^5$ no da 60°. L'accoppiamento misurato (5 mdeg, modo a 56.9 Hz) resta non modellato.
- **$T_{wy}$ ai modi in elevazione (diagnosi, 2026-11-02):** `analisi_lqg_twy_modi.py`. Il picco resta perché dipende dai poli dello stimatore: $q_f$ con $Q=V$ a $10^4$ muove solo il regolatore e fa da notch su $L$. A cuscinetti bloccati il tarato ha un polo a 57.8 Hz con $\zeta=0.001$ e amplifica 16 dB rispetto alla cascata, confermato su Dahl ad asse fermo. Con $V_f=10^6$ il picco scende di 6.4 dB senza costi visti; il controllore tarato non è stato cambiato.
- **Quando conviene l'LQG-PI (2026-11-02):** `analisi_lqg_benefici.py`, confronto con il P sull'encoder a guadagno alzato. Il vantaggio è il guadagno a parità di margine, grazie alla retroazione del polo del velocity loop: comandi a largo spettro e coppie sopra l'attrito a 3-7 Hz. Non c'è vantaggio su inseguimento regolare, coppie sotto l'attrito, rumore di misura. Non provati: lead, prefiltro, feedforward.
- **Vento (2026-11-02):** `analisi_lqg_benefici.py --vento`, raffiche di Davenport (livelli ipotizzati, A14). Rapporto cascata/LQG 6.0 in elevazione col tarato, 4.0 in azimut con $q_p=10^5$, 1.0 in azimut col tarato, perché il tuning ha fermato $q_p$ sulla corrente dei doppietti. Nel vento quel guadagno non chiede corrente in più. Da decidere: se tarare l'azimut su inseguimento e vento, con i comandi gestiti da limiti o preprocessore (A1). Attenzione: $q_p=10^5$ è instabile da 60° di elevazione senza gain scheduling.
- **Analisi per la tesi (2026-11-02):** `analisi_tesi_puntamento.py`, report §12 (`report_lqg_tuning/sec_thesis_pointing.tex`). Decisione dell'utente: i guadagni si fermano sui doppietti (A1 chiusa, D16); `config_drone_v2` si chiama `config_v3` nella tesi. Puntamento efficace con fascio da 8 cm, cascata → LQG-PI senza disturbo: circonferenza 13.8 → 29.1 %, dataset_7 config_v3 96.3 → 99.5 %, dataset_8 config_v3 53.5 → 61.0 %, EKF2 sotto il 3 %. Allineamento della stima su tutta la finestra; dai primi 5 s il limite di dataset_8 scende al 14 %.
- **Disturbo gaussiano con media 0.2 N m (2026-11-03, richiesta dell'utente):** `analisi_tesi_puntamento.py --media` (30 casi in più) e capitolo 11 (`analisi_lqg_circonferenza.py --senza-extra`, il vecchio "process noise" sostituito dal gaussiano corretto, vecchio riepilogo in `summary_process_noise_2026-11-03.json`). Le corse partono nell'equilibrio in cui l'integratore di velocita' del drive regge la media. Effetto: costante, -12.7 µrad in azimut e 0 in elevazione (flessione statica encoder-telescopio), +0.86 / +0.74 A; puntamento efficace al piu' -0.4 punti. Una coppia costante non e' il modo in cui contano i disturbi non modellati; le molle dei cavi identificate su V2 (0.93 / 1.66 N m/rad) non sono nel simulatore. §11.8 (benefici) usa ancora il vecchio process noise.
- **Simulazioni definitive per la tesi (2026-11-03):** `analisi_tesi_hil.py` -> `results_lqg_tuning/no_bearings/thesis_hil/`, figure in `Thesis-Writing/plots/hil/` (stile `style.py` della campagna), testo in `Thesis-Writing/sections/LQG_tests.tex` (Test Description, Results). Scenari: 0 nominale; 1 coppia gaussiana 0.3 N m rms 0.5-5 Hz su entrambi gli assi (3 realizzazioni); 2 Coulomb del plant x1.5; 3 combinato. Circonferenza e dataset7/8 con EKF2 e config_v3. Coppia: azimut uguale, polare 4.1-4.4 volte meglio con l'LQG-PI; dataset7 config_v3 eta 64.2 contro 89.9 %.
- **Traiettorie reali piu' veloci (2026-11-03):** `analisi_tesi_puntamento.py --velocita` -> `speed.json`, figura `speed`, report §12. Stesso percorso con config_v3 a 2, 3, 5, 10 volte la velocita' (0.19 -> 1.9 m/s), stimatore a 30 Hz reali, errore della stima della corsa originale. Distacco LQG-PI - cascata da 3-8 a 17-25 punti a 0.4-0.6 m/s, poi calano entrambi. L'elevazione rigida dell'LQG-PI risponde ai gradini della stima con picchi di corrente: oltre 18 A su dataset_8 da x3 e dataset_7 a x10 (il modello non satura). Servirebbe un comando piu' liscio delle stime tenute a 30 Hz.

**Punto da approfondire (utente):** cosa vuol dire nel dettaglio la linearizzazione armonica.
