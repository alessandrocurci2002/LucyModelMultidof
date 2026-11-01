# Troncamento modale alla Gawronski: elenco delle modifiche (2026-10-13)

> Tutti i numeri sono un **risultato da parametri segnaposto** (materiali, smorzamento,
> geometria polare 2-4, velocità di inseguimento): valgono per i dati attuali di
> `config.yaml`, non caratterizzano il sistema reale.

Documenti di riferimento: `SPEC_modello_gimbal_multidof.md`,
`MODELLO_unificato_gawronski_dedonatis.md` (§3-8, §10), `ANALISI_guyan_vs_gawronski.md`,
prototipo `tests/confronto_guyan_modale.py` (da cui vengono `zeri_da_modi` e la
flessibilità residua).

## 1. Modifiche per file

### Nuovi

| file | contenuto |
|---|---|
| `src/gimbal_mdof/modal.py` | Base modale del FEM COMPLETO di ogni asse con il payload congelato (`analizza_asse`, `modi_con_rigido`: modo rigido **1** non normalizzato, flessibili normalizzati in massa). Vettori d'aggancio c_e, c_m, c_L, c_P. Norme di Hankel (Gawronski 4.4-4.6). Scelta dei modi (taglio, Hankel, numero forzato, extra della validazione). Vettori statici residui (`base_modale`, opzione `residual_vectors`). Modello accoppiato **su base generica** `ModelloRitz` (eq. 5.1): basi completa (riferimento a 9 + 10 GdL), Guyan, Craig-Bampton, modale. Equilibrio e linearizzazione (eq. 7.2), stato a blocchi di Gawronski (2.25/2.32), confronto in banda (criteri T5), ciclo di validazione (`costruisci_sistema_modale`), controllo di maglia, flessibilità residua e zeri dell'asse isolato. M_p, C_p, g_p sono derivati con SymPy dal Lagrangiano del solo payload (`funzioni_payload`), con la cinematica di `coupled.py`, e compilati con `lambdify`. |
| `src/gimbal_mdof/modal_output.py` | CSV `modal_modes_<asse>.csv`, `modal_spectrum.csv`, `modal_vs_full.csv`; sezione di `report.md`; sezione LaTeX "Modal truncation (Gawronski)". |
| `tests/test_modal.py` | 23 test: T1-T9 più due controlli incrociati (base di Guyan = `coupled.py`; jacobiani del payload = forme chiuse di `Payload`). |
| `docs/CAMBIAMENTI_modale.md` | questo file. |

### Modificati

| file | modifica |
|---|---|
| `config.yaml` | Blocco `reduction.modal` aggiunto; `reduction.method: guyan → modal`; `frequency.f_max_Hz: 5000 → 2000` (commento sull'estensione automatica tolto). **Su tua richiesta:** il nodo 9 di azimut ora contiene il load "Tavola Rotante [DEF]" da Onshape (`frame: cad_cg`, `includes_polar_shafts: true`); il valore precedente, `J_kgm2: 1.0` segnaposto, è riportato in un commento. Nessun altro dato toccato. |
| `params.py` | `ParametriModali` (validazione di `reduction.modal`), `method: modal` predefinito. `InerziaNodale` con `frame: cad_cg` (J = I_zz,G + m d²) e `includes_polar_shafts`. Avviso `NODI_SENZA_INERZIA_<ASSE>`. Avviso Q7 declassato a nota quando il load comprende gli alberi polari. |
| `payload.py` | `J_ax = 0` (e `J_ax_nel_carico` per la tracciabilità) se il load di azimut comprende gli alberi polari. |
| `coupled.py` | Cinematica ed energie del payload estratte in `energie_payload` (stesse espressioni, riusate da `modal.py`); metodo `direzione_rigida`. |
| `simulate.py` | Le condizioni di moto rigido usano `acc.direzione_rigida` (in coordinate modali "tutte le ψ̇ = Ω₀" è "η̇_ψ0 = Ω₀"); `rtol`/`atol` come argomenti (default invariati 1e-10/1e-12). |
| `linearize.py` | `linearizza` e `costruisci_uscite` delegano a `ModelloRitz` quando serve. |
| `tf.py` | `esporta(..., extra=...)`: in forma modale salva anche Φ, ω, ζ, E_η, M̃, K̃, D̃ e lo stato di Gawronski (`A_gawronski`, `B_gawronski`, `C_gawronski`, `S_gawronski`); avviso segnaposto anche per `ModelloRitz`. |
| `analysis.py` | `sweep_theta0` e `frf_incrociate` accettano anche il sistema modale. |
| `plotting.py` | Griglia principale e secondaria, valori ed etichette con unità su ogni asse (richiesta sul PDF). Fasi del Bode allineate fra modelli (a un polo non smorzato `unwrap` sceglieva rami diversi, +360°). |
| `latex_report.py` | Tabelle 1-3 (larghezza e centratura). Nuova sezione 10 "Modal truncation (Gawronski)". Titolo, catena di modellazione, diagramma di flusso, tabella degli stadi e simulazioni aggiornati; Guyan resta come confronto. |
| `report.py` | Sezione modale; stato di Q4, Q7 e Q9 aggiornato; confronto di S1 con il modello a tutti i GdL. |
| `run_all.py` | Opzioni `--reduction modal\|guyan\|craig_bampton` e `--f-band`; griglie fino a `f_max` senza estensione automatica; pipeline modale; Bode `bode_modale_th*_om0`; esportazione in forma modale; S1 sul riferimento a tutti i GdL. |
| `README.md` | Riduzione modale, nuove tabelle, `node_inertias` con `frame: cad_cg`. Corretta la descrizione della terna 2, che era rimasta alla convenzione vecchia (Î = puntamento). |

`reduction.py`, `tests/config_riferimento.yaml` e i 92 test preesistenti **non** sono
stati modificati.

## 2. Test di accettazione

| richiesta | test | esito |
|---|---|---|
| T1 | `test_T1_base_modale` (2 assi × 2 θ_ref) | K1 = 0, ΦᵀM̄Φ = diag(J_tot, 1, …) e ΦᵀKΦ = diag(0, ω²) a 1e-10, prima colonna 1, c_e0 = c_m0 = c_agg,0 = 1, ΦᵀβKΦ = diag(2ζω) |
| T2 | `test_T2_residuo_non_lineare` (θ_ref = 0 e 25°) | residuo modale = Φᵀ · residuo fisico entro 1e-10 su 20 stati casuali; anche l'energia coincide |
| T3 | `test_T3_frf_tutti_i_modi`, `test_T3_blocchi_di_gawronski`, `test_T3_jacobiano_numerico` | FRF modale e stato di Gawronski = FRF fisica entro 1e-8 (1 Hz - 2 kHz, tutti i punti di lavoro); poli lenti coincidenti; blocchi (2.25)/(2.32) esatti nel caso disaccoppiato; la (7.2) coincide con lo jacobiano numerico |
| T4 | `test_T4_limite_rigido` (3 punti) | solo i modi rigidi = modello a 2 GdL linearizzato entro 1e-10 |
| T5 | `test_T5_troncato_in_banda`, `test_T5_ciclo_di_validazione` | criteri rispettati in tutti i punti dopo il ciclo |
| T6 | `test_T6_vettori_residui` | zeri dei canali diretti entro 1e-3, frequenze dei modi tenuti spostate < 1e-3 (in pratica ~1e-11) |
| T7 | `test_T7_frequenze_esatte` | frequenze dei modi tenuti = FEM completo entro 1e-10 |
| T8 | `test_T8_ordine_dei_modi` | rigidi per primi (azimut, polare), poi solo flessibili crescenti; numero di stati coerente |
| T9 | `test_T9_S1_modale_contro_lineare`, `test_T9_S3_conservazione`, `test_T9_S1_modale_contro_completo` | S1 < 1 % sui canali di azimut; S3: E e p_ψ conservati sotto 1e-8; S1 modale contro tutti i GdL entro 1e-3 del picco |

Più due controlli incrociati: `test_base_guyan_uguale_a_coupled` (la base di Guyan nel
modello generico riproduce residuo, energia e M₀, C₀, K₀ del modello SymPy di
`coupled.py` entro 1e-12) e `test_payload_linearizzato` (γ, k_c, k_g, M₁₂).

**pytest:** 115 test superati (92 preesistenti invariati + 23 nuovi), ~55 s.

### Dove ho dovuto interpretare la richiesta (da rivedere)

1. **T3, griglia da 1 Hz.** Sotto ~1 Hz la FRF del modello FISICO in doppia precisione
   ha un pavimento d'errore ε|K|/(ω²J) ≈ 1e-9…3e-8: voci di K ~1e7 N·m/rad contro ω²J ~ 1,
   con il polo a pendolo polare non smorzato a 0.16 Hz. L'ho verificato con mpmath a 40
   cifre: a 0.164 Hz numpy sbaglia di 3.3e-8 sul modello fisico e di 7e-9 sul modale.
   La FRF si confronta quindi da 1 Hz a 2 kHz, escludendo ±5 % attorno ai poli non
   smorzati. La dinamica lenta si verifica confrontando direttamente i poli sotto 10 Hz.
   Il doppio polo nullo di azimut del modello fisico vale ±6e-5 j rad/s, per la "molla a
   terra" nata dagli arrotondamenti di K·1: da qui `atol` = 1e-4 rad/s.
2. **T6 sul modello predefinito** (modi della validazione + vettori residui). Con il
   solo taglio, sulla configurazione di prova, uno zero *dinamico* polare (938 Hz) resta
   a 1.3e-3: dipende dall'inerzia dei modi scartati, non dalla loro flessibilità.
3. **T9 contro tutti i GdL**: il riferimento si integra con rtol 1e-8 e atol 1e-10. Il FEM
   completo ha modi fino a ~260 kHz (configurazione di prova) e il doppietto discontinuo
   li eccita davvero, quindi DOP853 è limitato dalla stabilità. Con 1e-10 serviva circa
   1 s di calcolo per ogni ms simulato; anche Radau e BDF sono lentissimi alla stessa
   tolleranza.
4. **`n_flex` forzato** prende i primi n modi in ordine di frequenza. La scelta dei modi
   si fa nel `theta_ref_deg` di config e resta la stessa, per indice, a ogni θ₀; le
   frequenze cambiano di < 0.01 % in ±45°.

## 3. Risultati con `config.yaml` attuale (consistente, θ_ref = 0)

| | azimut | polare |
|---|---|---|
| J_tot [kg·m²] | 2.1256 | 0.88324 |
| modi flessibili tenuti [Hz] | 752.2 (taglio), **4595.1 (validazione)** | 335.6, 3557.2 (taglio) |
| primo scartato [Hz] | 21 344 | 14 568 |
| flessibilità residua collocata [rad/(N·m)] | 7.1e-8 (9e-4 del totale) | 6.7e-8 (7e-4) |

* **Stati: 12** contro i 14 di Guyan. La stima a priori della richiesta era ~10 stati
  (1 modo azimutale e 2 polari). Il ciclo di validazione ha aggiunto il modo azimutale 2
  (4595 Hz, sopra il taglio di 4 kHz). Senza di esso lo zero carico-molla di
  τψ→ψe (12.49 Hz) si sposta dell'1.4 %, e T5 (zeri entro 1 %) non passa: quel modo porta
  gran parte della flessibilità statica fra motore e forcella, perché i nodi motore ed
  encoder sono quasi privi di massa (vedi §4.a). Con `static_correction:
  residual_vectors` lo stesso zero torna esatto anche senza il modo 2.
* **T5** in tutti i 10 punti di lavoro: canali diretti ≤ 0.018 dB e 0.061°; zeri entro
  9e-5; matrice 2×2 entro 2.2e-4. **Guyan**, con gli stessi criteri: 22 dB / 157°
  (risonanza azimutale a 834 Hz invece di 752 Hz), zeri sbagliati fino al 14 % (polare
  504.5 → 577 Hz).
* **Zeri dell'asse isolato** (puro → con correzione statica): azimut τ→motore 12.4968 →
  12.491 Hz (completo 12.491); polare τ→motore 17.1847 → 17.1788 (completo 17.1788).
* **Controllo di maglia** (concentrata contro consistente): −0.02 %, −0.48 % (azimut);
  −0.01 %, −0.79 % (polare). Nessun avviso.
* **Simulazioni** (modello modale, θ₀ = 0°, cioè il punto con k_g > 0 più vicino a +15°):
  - S1 contro il lineare modale: 9.9e-7 su ψ_e, 1.2e-5 su ψ_m, 3.8e-10 su ψ_L. I canali
    polari valgono ~0.9, ma come già documentato rispondono solo per i prodotti d'inerzia,
    quindi non misurano la linearizzazione.
  - S2: 5.8e-7 su ψ_e.
  - S3: derive di E e p_ψ pari a 2e-14 e 2.6e-12.
  - S1 sul modello completo, nei primi 0.1 s: scarto massimo, rispetto al picco, 1.6e-5 su
    ψ_e, 3.3e-4 su ψ_m e 6e-9 su ψ_L. Lo scarto su ψ_m viene dai modi sopra il taglio, che
    il doppietto discontinuo eccita sul nodo motore, leggero.

## 4. Incoerenze e osservazioni (segnalate, non corrette in silenzio)

a. **Dati.** `polar.node_inertias` è tutto commentato, e nell'azimut mancano l'inerzia di
   rotore al nodo 5 e di encoder al nodo 1. I nodi motore ed encoder hanno quindi solo la
   massa degli alberi. Ora c'è un avviso automatico (`NODI_SENZA_INERZIA_AZIMUT`,
   `NODI_SENZA_INERZIA_POLARE`). È la causa del modo azimutale aggiunto dalla validazione.
b. **Percorsi.** La richiesta indica `docs/MODELLO_unificato…`, `docs/ANALISI…` e
   `scripts/confronto_guyan_modale.py`: i due documenti sono nella radice del progetto,
   il prototipo in `tests/`. Non li ho spostati.
c. **Numero di test.** La richiesta dice 91 test esistenti: erano 92.
d. **`MODELLO_unificato` §5 e §7.1** usano la vecchia convenzione su θ (dall'orizzonte,
   positivo verso il basso): g_p = −m g (d_x cos θ + d_z sin θ), k_g = m g (d_x sin θ −
   d_z cos θ). Il codice usa θ > 0 che alza il puntamento, con k_g = −m g (d_x sin θ +
   d_z cos θ). Ho ricavato M_p, C_p e g_p con SymPy dalla cinematica di `coupled.py`, e
   il test `test_payload_linearizzato` li confronta con `Payload`.
e. **`MODELLO_unificato` §3-4** costruisce Φ dalle catene già ridotte con Guyan (come
   segnala `ANALISI` §5.1): qui Φ viene dal FEM completo, come richiesto.
f. **`Sistema.isolato_completo`** mette nel polare anche k_g(θ₀), mentre i modi vanno
   calcolati con k_g = 0. `modal.py` usa la massa di riferimento senza k_g, e k_c + k_g
   entrano come accoppiamento di rango 1 (variante predefinita, §7.2).
g. **`reduction.method: craig_bampton`** prima non era utilizzabile nel modello
   accoppiato SymPy (le matrici di Craig-Bampton non sono tridiagonali); ora passa per
   `ModelloRitz`.
h. **Smorzamento C = βK**: i modi sopra ~20 kHz hanno ζ > 0.3, fino a ζ > 1 (azimut modo
   6: 1.11; polare modo 4: 1.46). È una conseguenza del modello proporzionale alla
   rigidezza: non tocca la banda, ma va tenuto presente se si alza il taglio.
i. **Correzione statica e zeri dinamici.** I vettori residui rendono esatti gli zeri
   "statici", ma possono spostare leggermente quelli dinamici. Sulla configurazione di
   prova lo zero polare a 938 Hz passa da 9.5e-6 a 2.7e-4: aggiungono flessibilità senza
   la relativa inerzia.
j. **Forcella (Onshape).** Hai chiesto se fosse corretto che il load contenga gli alberi
   polari. Non lo era del tutto: il modello aggiungeva già il loro momento attorno
   all'asse di azimut (J^ax) al nodo L, quindi contava gli alberi due volte. Il
   nodo 9 e J^ax stanno sullo stesso nodo, quindi sottrarre gli alberi dal load e poi
   riaggiungerli come J^ax equivale a non aggiungere J^ax. Il totale diventa esattamente
   il valore CAD, che include anche gli Steiner m y_c² veri: prima mancavano, perché Q7
   (`y_node1_from_o2_mm`) è nullo. J^ax calcolato vale 1.07e-3 kg·m², non aggiunto.
   L'inerzia degli alberi attorno all'asse polare resta nel FEM polare e non è contata due
   volte. Da verificare:
   - che "Tavola Rotante [DEF]" non contenga anche parti dell'albero di azimut (elementi
     1-8), altrimenti anche quelle sarebbero contate due volte;
   - Onshape mostra L_yy e L_zz con sole 4 cifre (1.592e6 e 1.248e6 kg·mm²): conviene
     aumentare la precisione di visualizzazione;
   - il baricentro della forcella è a 12.2 mm dall'asse di azimut: nel modello torsionale
     entra solo m d² = 0.0042 kg·m² (J passa da 1.248 a 1.2522 kg·m²).
k. **README**: la descrizione della terna 2 era rimasta alla convenzione vecchia
   (Î = puntamento); l'ho corretta (K̂ = puntamento, Î verso l'alto).
l. **`theta_ref`**: `MODELLO_unificato` §12 suggerisce −45° per le simulazioni non
   lineari; la richiesta fissa 0° (centro del campo operativo). Uso 0°, dalla
   configurazione. La (5.1) è esatta per qualunque scelta.

## 5. Esito di `run_all.py`

* `run_all.py --quick`: completato senza errori in 3 min 43 s; PDF compilato.
* `run_all.py` (completo, ultimo eseguito, quindi sono suoi i file in `results/`):
  completato in 4 min 52 s. Log: 12 stati, 1 modo aggiunto dalla validazione, T5
  rispettato. PDF di 30 pagine, 0 errori LaTeX, nessun riferimento irrisolto.
* Il tempo cresce rispetto a prima (~2 min) soprattutto per S1 sul modello a tutti i GdL
  (finestra di 0.1 s) e per i Bode modale/completo/Guyan.
* Restano 5 righe fuori margine nel PDF (da 5 a 47 pt), tutte già presenti prima di
  questa modifica (paragrafi con formule lunghe nelle sezioni 2, 7 e 11).

## 6. Aggiornamento del 2026-10-15: PDF senza Guyan e terne

* **`latex.include_guyan`** (nuovo, `config.yaml`, ora `false`; se manca vale `true`, quindi
  la configurazione congelata dei test non cambia). Con `false` il PDF non contiene:
  - la sezione di Guyan e la catena a tre inerzie in forma chiusa;
  - il modello accoppiato SymPy a 7 GdL, le FdT simboliche delle catene e le figure
    completo/ridotto;
  - le colonne e la curva "Guyan" nella sezione modale.

  Al loro posto ci sono:
  - la sezione *Coupled model* scritta con tutti i GdL del FEM (eq. `eq:full_eom`, termini
    del telescopio M_p, C_p q̇_p, g_p);
  - la sezione *Operating point and linearisation* (γ, k_c, k_g e la tabella dei punti di
    lavoro);
  - i risultati del modello lineare modale: errore di misura, frequenze in funzione di θ₀,
    modi del FEM completo, poli e zeri.

  `report.md` e i CSV continuano a riportare il confronto con Guyan.
* **Figure**:
  - `errore_misura_*` ora usa il modello principale (modale), non più Guyan;
  - nuova figura `modi_modale_<asse>`;
  - il Bode modale contiene la curva di Guyan solo con `include_guyan: true`.
* **Terne**: aggiunto il paragrafo "Why Ĵ = −ĵ" nella sezione *Reference frames*. Codice e
  convenzione non cambiano.
* **Corretto un difetto preesistente**: `latex_report.L()` non convertiva i simboli delle
  matrici SymPy (le matrici mutabili non sono `sp.Basic`). Per questo nel PDF la matrice
  M(q) usciva con I_xx, I_zz invece di I_KK, I_II.
* **Esito**: pytest 115 test superati; `run_all.py` completato senza errori in 4 min 55 s;
  PDF di 26 pagine, 0 errori LaTeX, nessun riferimento irrisolto. Restano 2 righe
  leggermente fuori margine (14 e 37 pt), preesistenti.

## 7. Aggiornamento del 2026-10-15: modello a parametri concentrati (3 + 4 GdL)

Richiesta: azimut 1-5-9, polare 1-P-7-11; elementi equivalenti con inerzia somma e rigidezza
in serie (1/k = Σ 1/k_i); stesso accoppiamento dinamico fra gli assi del modello accoppiato.

* **Nuovi moduli**
  - `concentrato.py`: elementi equivalenti e masse concentrate ai nodi (metà di ogni
    elemento adiacente più l'inerzia nodale; le inerzie dei nodi intermedi, oggi nessuna,
    vanno nell'elemento). Il sistema è costruito con la catena di `coupled.py`, quindi
    Lagrange, linearizzazione, FdT e simulazioni sono gli stessi. Le **equazioni di
    Newton-Eulero** sono ricavate con SymPy in modo indipendente (bilancio dei momenti su
    ogni nodo; equazioni di Eulero del telescopio attorno al punto fisso o_2) e verificate
    identiche a Lagrange.
  - `concentrato_output.py`: sezione LaTeX con le equazioni per GdL, in forma compatta e
    in forma esplicita con le funzioni trigonometriche, e verifica della forma compatta.
  - `concentrato_analisi.py`:
    - FdT simboliche delle catene, con le antirisonanze in forma chiusa: √(k_ψ,2/J_ψ,L),
      √(k_ψ,1/J_ψ,e), √(k_θ,3/J_θ,r);
    - linearizzazione in tutti i punti di lavoro e criteri T5 contro il modello completo.
  - `simulatore_concentrato.py`:
    - simulatore non lineare M(θ_p) q̈ = E τ − h(q, q̇), dalle equazioni di
      Newton-Eulero, con la stessa interfaccia di `simulate.py`;
    - scenario S4 a grande ampiezza;
    - esportazione del codice autonomo: `results/simulatore_concentrato/ne_concentrato.py`
      (verificato) e `ne_concentrato.m` più `ne_concentrato_parametri.m` (generati dalle
      stesse espressioni, **non eseguiti**: MATLAB/Octave non disponibili qui).
  - `tests/test_concentrato.py`: 8 test.
* **Modifiche**
  - `run_all.py`: blocco del modello concentrato (FdT, Bode, sensibilità a θ₀,
    S1/S3/S4a-c, CSV, esportazioni).
  - `latex_report.py`: sezione `sec:lumped`.
  - `plotting.py`: `storie_grande_ampiezza`.
  - `report.py`: campo `extra_md`.
* **Risultati** (dati attuali, segnaposto):
  - **Frequenze degli assi isolati**, concentrato contro FEM completo: azimut 774 contro
    752 Hz e 1847 contro 4595 Hz; polare 379 contro 336, 943 contro 3557 e 3343 contro
    14568 Hz.
  - **Cause degli scarti**: la divisione a metà dell'inerzia porta metà dei dischi pesanti
    (polare elementi 2 e 9) su nodi lontani, e mancano le masse intermedie.
  - **Antirisonanza carico-molla** (12.49 Hz a θ₀ = 0): esatta.
  - **Criteri T5 non rispettati** sui modi strutturali (fino a 27 dB e 173° vicino alle
    risonanze spostate): è il limite atteso del modello a 3 + 4 nodi.
  - **S4**: lo scarto fra non lineare e lineare cresce con l'ampiezza (θ₀ = −45°: 1.2 % a 2°,
    8.5 % a 10°). A θ₀ = 0° un'oscillazione di 10° supera il limite di stabilità del
    telescopio sbilanciato (instabile sopra ~+4°) e il telescopio arriva a ~40°, mentre il
    modello lineare prevede un'oscillazione limitata.

## 8. Aggiornamento del 2026-10-16: confronto con le misure sperimentali

* **Dati**: `sim/LucyIdentification/*.csv` (Kollmorgen WorkBench, misura "Plant"); la
  cartella `Validation` è esclusa (il parser di `config.yaml` la rifiuta).
  `Elevation_new_BodeMeasurement.csv` ha gli stessi dati di misura di
  `Elevation_BodeMeasurement.csv` (cambia solo la taratura salvata) e non è ripetuto.
* **Unità**: la pianta di WorkBench è velocità di retroazione / corrente iniettata, in
  (rad/s)/A. Lo conferma il tratto rigido dell'elevazione a 10 Hz: K_t/(Jω) del modello dà
  −46.6 dB, la misura −46.8 dB. Il modello è convertito come G = jω K_t θ/τ, con
  K_t = 0.258 N m/A (TBM(S)-7646-A a 25 °C, A rms).
* **Modello**: solo il 7 GdL concentrato, linearizzato all'elevazione della prova con
  Ω₀ = 0, ai nodi encoder e motore. Per la prova senza telescopio il payload è tolto.
* **Nuovi moduli e file**:
  - `sperimentale.py`, `sperimentale_output.py`;
  - `tests/test_sperimentale.py` (4 test, con CSV sintetici);
  - sezione `experimental` in `config.yaml` e sezione `sec:experimental` del PDF;
  - figure `exp_*.pdf`, tabella `experimental_comparison.csv`.
* **`plot_bode.py`** (in `LucyIdentification`):
  - stile delle figure del report: palette, griglia principale e secondaria, unità, PDF e
    PNG a 300 dpi;
  - asse da 0.1 Hz, con una fascia grigia dove non ci sono dati;
  - opzione `--all`;
  - rigenerati i `*_bode.png/.pdf` della cartella.
* **Ipotesi da confermare**:
  - elevazione della prova di elevazione assunta θ₀ = 0° (non registrata);
  - "Po_30" interpretato come elevazione +30°;
  - corrente in A rms.
* **Risultati principali** (dettaglio nel PDF, §11):
  1. Fra 5 e 300 Hz il modello riproduce l'azimut: tratto rigido, antirisonanza
     carico-molla (11.9 Hz misurati contro 12.5 Hz) e il suo spostamento con l'elevazione.
  2. Sotto 5 Hz c'è una rigidezza verso terra che il modello non ha (risonanza a 3.8 Hz,
     k₀ ≈ 1.2e3 N m/rad). Dipende dall'ampiezza (4.2 Hz con 2.36 A): è l'attrito dei
     cuscinetti in regime di pre-strisciamento.
  3. Primo modo strutturale dell'azimut: 490 Hz misurati contro 774 Hz. Con le rigidezze del
     modello servirebbe un'inerzia lato motore di 9.9e-4 kg m² invece di 4.0e-4: mancano il
     rotore (8.19e-5) e il disco encoder, e probabilmente c'è cedevolezza nelle giunzioni.
  4. Livello −4.9 dB fra 15 e 300 Hz: da controllare prima la scala di coppia (A picco
     contro A rms: −3 dB; avvolgimento caldo: −0.8 dB).
  5. Coppia risonanza/antirisonanza a 37/45 Hz in azimut senza corrispondenza nel modello
     (forcella supposta rigida, base, cavi).
  6. In elevazione l'antirisonanza a 19 Hz corrisponde al canale coppia→motore (17.2 Hz),
     non a quello coppia→encoder; fra 50 e 250 Hz ci sono modi non modellati.
  7. Senza telescopio l'asse si comporta come una molla (k₀ ≈ 2.8e3 N m/rad) con
     J ≈ 0.036 kg m², contro 1.8e-3 del modello.
