# Modello dinamico multi-GdL del pointer LUCY

Modello principle-based del pointer alt-azimutale LUCY, con la sua validazione sulle misure. La catena di
modellazione è questa:

1. FEM torsionale degli alberi di azimut e polare;
2. riduzione per troncamento modale alla Gawronski (Guyan e Craig–Bampton come confronto);
3. accoppiamento dei due assi tramite l'energia cinetica;
4. modelli concentrati da 7 a 9 GdL;
5. stima modale ai minimi quadrati sulle FRF misurate;
6. attrito dei cuscinetti di Dahl, identificato su frequenza e tempo;
7. verifica in anello chiuso con il modello del drive AKD2G.

Implementa `SPEC_modello_gimbal_multidof.md`. È il secondo modello della catena dopo il modello rigido
`../gimbal_2dof`, che viene **riusato** (dati CAD del telescopio, preprocessing inerziale, convenzioni, stile delle
figure). Il modello in anello chiuso e gli studi di controllo stanno in [`closedLoopModel/`](closedLoopModel/README.md).

> **Dati segnaposto.** `config.yaml` contiene dati geometrici, materiali, inerzie e smorzamenti segnaposto; anche H7
> ha ancora 18 voci `placeholder: true`. Gli output che ne dipendono sono marcati "risultato da parametri segnaposto".

## Indice

- [Stato del progetto](#stato-del-progetto)
- [Installazione](#installazione)
- [Avvio rapido](#avvio-rapido)
- [Dati esterni e dipendenze](#dati-esterni-e-dipendenze)
- [Struttura della cartella](#struttura-della-cartella)
- [Versioni del modello](#versioni-del-modello)
- [Script di analisi](#script-di-analisi)
- [Report](#report)
- [Risultati principali](#risultati-principali)
- [Risultati di `run_all.py`](#risultati-di-run_allpy)
- [Come inserire la geometria](#come-inserire-la-geometria)
- [Codice sorgente](#codice-sorgente)
- [Test](#test)
- [Versionamento](#versionamento)
- [Limiti noti](#limiti-noti)
- [Decisioni e scostamenti dalla specifica](#decisioni-e-scostamenti-dalla-specifica)

## Stato del progetto

Aggiornato al 2026-11-01.

| fase | stato | dove |
|---|---|---|
| FEM torsionale e riduzione modale | fatto | `src/gimbal_mdof/`, `run_all.py` |
| Modelli concentrati H1–H7 | fatto; H7 è il riferimento | `config_hub*.yaml` |
| Confronto con le FRF misurate | fatto | `development_report/`, `analisi_20261028.py` |
| Modi della forcella, stima modale | fatto | `analisi_forcella_modale.py`, `report_fork_modes/` |
| Cuscinetti di Dahl, identificazione ibrida | fatto | `analisi_cuscinetti_dahl_tempo.py`, `report_dahl_bearings/` |
| Anello chiuso con il modello ibrido | fatto, due fenomeni non riprodotti | `analisi_anello_chiuso_ibrido.py` |
| Osservatore di disturbo con giroscopio | simulato | `closedLoopModel/report_dob/` |
| LQG alla Gawronski | **solo audit di progetto, non implementato** | `closedLoopModel/LQG_GAWRONSKI_AUDIT.md` |

## Installazione

Servono Python ≥ 3.10 e `latexmk` con `pdflatex` per i report. La cartella `../gimbal_2dof` deve esistere accanto a
questa.

```bash
cd sim/gimbal_multidof
python3 -m venv .venv
.venv/bin/pip install -e ../gimbal_2dof     # modello a 2 GdL, editable
.venv/bin/pip install -e ".[dev]"
```

## Avvio rapido

```bash
.venv/bin/python -m pytest                                 # test
.venv/bin/python run_all.py --quick --no-latex             # catena FEM su griglie ridotte, senza LaTeX
.venv/bin/python run_all.py --config config_hub_telescopeonring_alhub_fit_ground_aztel.yaml   # catena su H7
.venv/bin/python -I analisi_forcella_modale.py             # stima modale ultima (richiede ../LucyIdentification)
```

Opzioni di `run_all.py`:

| opzione | effetto |
|---|---|
| `--config FILE` | configurazione; i risultati vanno nella sua `results_dir` |
| `--hub` | scorciatoia per `config_hub.yaml` |
| `--quick` | griglie ridotte e simulazioni brevi (circa 4 min) |
| `--no-latex` | non genera né compila il LaTeX |
| `--theta0 ...`, `--omega0 ...` | punti di lavoro [deg] e velocità di inseguimento [deg/s] |
| `--mass-model lumped\|consistent\|all` | modello di massa |
| `--mass-form guyan\|guyan_lumped\|all` | forma di massa ridotta |
| `--reduction modal\|guyan\|craig_bampton` | riduzione del modello principale |
| `--f-band HZ` | banda di interesse |

## Dati esterni e dipendenze

Questo repository non è autosufficiente. Le misure e il modello a 2 GdL stanno fuori, a percorsi relativi:

| percorso | contenuto | necessario per |
|---|---|---|
| `../gimbal_2dof` | pacchetto `gimbal`, modello rigido, repository git separato | installazione, test, `run_all.py` |
| `../LucyIdentification` | misure WorkBench, sequenze, file `.akd`, datasheet; **non è sotto git** | tutti gli `analisi_*.py`, `closedLoopModel`, confronto sperimentale |
| `../kollmorgen_control_report` | report sulle maglie di controllo Kollmorgen | solo riferimento per la struttura del drive |

Sottocartelle di `../LucyIdentification` usate qui:

| cartella | contenuto |
|---|---|
| `VAntiRes/` | prova con 0–12 viti anello-telescopio e giochi radiali (19/10/2026) |
| `VClosedLoop/` | Bode, sequenze a step, `.akd` di ottobre |
| `VClosedLoop/V2/` | doppietti da 0.1 a 15°, un asse per file |
| `VClosedLoop/VTesi/` | FRF con 14 viti (23/10/2026) |
| `Poli_Polar_identification/` | misure polari delle varianti di montaggio |
| `Validation/` | sequenze di settembre, usate solo da `closedLoopModel` |

Esistono anche dipendenze **interne** fra cartelle di risultati: alcuni script leggono i `summary.json` di passi
precedenti. Sono indicate nella sezione [Script di analisi](#script-di-analisi); prima di cancellare o ignorare una
cartella `results_*`, verificare che nessuno la legga.

## Struttura della cartella

```
gimbal_multidof/
├── README.md, pyproject.toml, requirements.txt, .gitignore
├── SPEC_modello_gimbal_multidof.md              specifica del modello e delle domande Q1, Q2, ...
├── MODELLO_unificato_gawronski_dedonatis.md     rigido + flessibile: De Donatis come proiezione sui modi rigidi
├── ANALISI_guyan_vs_gawronski.md                confronto fra le riduzioni
├── docs/CAMBIAMENTI_modale.md                   modifiche del passaggio al troncamento modale (2026-10-13)
├── config.yaml, config_hub*.yaml                configurazioni: una per versione del modello
├── run_all.py                                   catena FEM completa, una configurazione alla volta
├── analisi_*.py                                 analisi puntuali sulle misure
├── src/gimbal_mdof/                             libreria
├── tests/                                       test e configurazione congelata
├── results/, results_hub*/                      uscite di run_all.py, una cartella per configurazione
├── results_fork_modes*/, results_20261028/, results_V2_coupling/    uscite di analisi_*.py
├── saved_fork_modes_modal_20261031_cad20_kt10/  copia di sicurezza; è un ingresso di altri script
├── development_report/, report_fork_modes/, report_dahl_bearings/   report LaTeX generati
└── closedLoopModel/                             anello chiuso, identificazione, controllo
```

## Versioni del modello

Ogni versione è una copia della precedente con modifiche minime, descritte nell'intestazione del file. Per
tracciabilità, ognuna scrive in una propria cartella di risultati.

| versione | configurazione | `results_dir` | modifica |
|---|---|---|---|
| base (7 GdL) | `config.yaml` | `results/` | telescopio rigido in P |
| H1 | `config_hub.yaml` | `results_hub/` | mozzo polare e attacco cedevole k_b del telescopio |
| H2 | `config_hub_telescopeonring.yaml` | `results_hub_telescopeonring/` | tratto del tubo nell'anello solidale al mozzo |
| H3 | `config_hub_telescopeonring_alhub.yaml` | `results_hub_telescopeonring_alhub/` | mozzi laterali in alluminio (errore di materiale nel CAD) |
| H4 (8 GdL) | `config_hub_telescopeonring_alhub_stiffring.yaml` | `results_hub_telescopeonring_alhub_stiffring/` | anello irrigidito |
| H5 | `config_hub_telescopeonring_alhub_fit.yaml` | `results_hub_telescopeonring_alhub_fit/` | due parametri liberi ritarati sulla prova delle viti |
| H6 | `config_hub_telescopeonring_alhub_fit_ground.yaml` | `results_hub_telescopeonring_alhub_fit_ground/` | vincolo a terra k(1+jη), k_b = 1.16e5 N m/rad |
| **H7 (9 GdL)** | `config_hub_telescopeonring_alhub_fit_ground_aztel.yaml` | `results_hub_telescopeonring_alhub_fit_ground_aztel/` | telescopio appeso al load di azimut |

H7 è la base predefinita di tutti gli script recenti. In `closedLoopModel` gli studi di settembre e ottobre usano H6.

## Script di analisi

Si eseguono dalla radice con `.venv/bin/python -I <script>`. Salvo diversa indicazione, la configurazione predefinita
è H7. In ordine cronologico:

| script | domanda | uscita | opzioni utili |
|---|---|---|---|
| `analisi_nodo_azimut.py` | Il telescopio appeso al load spiega i ~50 Hz in azimut? | `results_hub_..._aztel/` | `--config`, `--base`, `--k-sweep` |
| `analisi_rigidezza_assi.py` | Il giunto ha la stessa rigidezza in polare e in yaw? | `results_hub_..._aztel/` | `--config`, `--k` |
| `analisi_20261028.py` | Confronto delle misure del 28/10 con H7, rigidezza k_b | `results_20261028/` | `--config`, `--k-sweep` |
| `analisi_v2_accoppiamento.py` | Accoppiamento fra gli assi e bilanciamento dai doppietti | `results_V2_coupling/` | |
| `analisi_modi_forcella.py` | Modi della forcella eccitati dallo statore e letti dall'encoder | `results_fork_modes/` | `--f-min`, `--f-max`, `--modi`, `--cartella` |
| `analisi_attrito_cuscinetti.py` | Dahl sui cuscinetti spiega la risonanza a 6–9 Hz? | `results_fork_modes_2masses/bearing_friction/` | |
| `analisi_forcella_modale.py` | Stima modale su Re/Im con vincoli fisici | `results_fork_modes_modal/` | `--zeta-max`, `--cad`, `--cad-J`, `--kt`, `--cartella`, `--partenza` |
| `analisi_forcella_lugre.py` | Stessa stima con LuGre o Dahl al posto di k_g(1+jη) | `results_fork_modes_lugre/`, `results_fork_modes_dahl/` | `--attrito lugre\|dahl` |
| `analisi_cuscinetti_dahl_tempo.py` | FRF con Dahl dal tempo, stima ibrida di T_c e σ₀ | `results_fork_modes_dahl/time_domain_bearings/` | |
| `analisi_tempo_dahl_ibrido.py` | Modello ibrido nel tempo, ingresso corrente misurata | `results_fork_modes_dahl/time_domain_check/` | |
| `analisi_anello_chiuso_ibrido.py` | Modello ibrido in anello chiuso, ingresso riferimento di posizione | `results_fork_modes_dahl/closed_loop_check/{stima,akd}/` | `--guadagni stima\|akd` |

### Ordine di esecuzione della catena attuale

```
closedLoopModel/v2_identification.py        -> closedLoopModel/results_v2_identification/summary.json
analisi_modi_forcella.py --cartella results_fork_modes_2masses (*)
analisi_forcella_modale.py                  (parte da results_fork_modes_modal_cad20/, confronta con _2masses/)
saved_fork_modes_modal_20261031_cad20_kt10/ (copia manuale della stima approvata)
analisi_forcella_lugre.py --attrito dahl    (legge saved_... e results_v2_identification)
analisi_cuscinetti_dahl_tempo.py            (legge results_fork_modes_dahl e results_v2_identification)
analisi_tempo_dahl_ibrido.py                (legge time_domain_bearings/summary.json)
analisi_anello_chiuso_ibrido.py             (idem, più lucy_cl e VClosedLoop)
report_dahl_bearings/build_dahl_report.py
```

(\*) Il comando che ha prodotto `results_fork_modes_2masses/` non è registrato; il `summary.md` è quello di
`analisi_modi_forcella.py`. Anche altre cartelle hanno una provenienza parziale:

| cartella | provenienza |
|---|---|
| `results_fork_modes_modal_zeta001`, `_zeta005_J15`, `_zeta005_J20`, `_cad20` | passi precedenti di `analisi_forcella_modale.py`, con `--cartella` e limiti diversi |
| `results_fork_modes_modal_diag_*`, `_zeta001_diag_cad50` | diagnostiche con limiti allargati; solo `_diag_cad50` ha il comando documentato |
| `results_fork_modes_time_domain`, `results_fork_modes_k2_test` | prove una tantum, **senza script nel repository** |

## Report

Testo in inglese, figure e numeri generati dai risultati. I tre script di report leggono le uscite esistenti e non
rifanno le stime, tranne `build_development_report.py`, che ricalcola dai file di configurazione.

| report | pagine | contenuto | rigenerazione, dalla radice |
|---|---|---|---|
| `development_report/development_report.pdf` | 31 | da 7 a 8 GdL, prova delle viti, H5 e H6 | `.venv/bin/python development_report/build_development_report.py` |
| `report_fork_modes/fork_modes_report.pdf` | 15 | stima modale dei modi della forcella | `cd report_fork_modes && ../.venv/bin/python -I build_fork_modes_report.py` |
| `report_dahl_bearings/dahl_bearings_report.pdf` | 14 | Dahl, stima ibrida, verifica in frequenza, tempo e anello chiuso | `cd report_dahl_bearings && ../.venv/bin/python -I build_dahl_report.py` |
| `results*/latex/multidof_fem_model.pdf` | | analisi FEM con la geometria in TikZ | `run_all.py` senza `--no-latex` |

`--no-latex` è disponibile su tutti e tre gli script. `report_fork_modes/previous_logfit/` conserva la versione
precedente del report (stima con errore logaritmico).

## Risultati principali

Sintesi dei report; i dettagli e i limiti sono nei PDF.

**Il modello ha guidato la diagnosi dell'assemblato.**

- Prova con 0–12 viti anello-telescopio: da 4 a 8 viti k_b cresce di circa 9 volte e l'antirisonanza passa da
  17.5 Hz a circa 58 Hz; giochi non uniformi la spostano del 5–15 %.
- Le risonanze a ~56 e ~99 Hz sono modi della struttura di forcella e statore, letti dall'encoder. Un residuo modale
  negativo produce zeri a fase non minima.
- In azimut, sopra il modo della forcella, il livello misurato chiede k₂ circa 1.6 volte il CAD. L'82 % della
  cedevolezza di k₂ viene dal tratto Ø16 × 30.88 mm, da verificare sull'hardware.

**Cuscinetti di Dahl, stima ibrida** (T_c dai doppietti V2, σ₀ dalle FRF):

| asse | σ₀ [N m/rad] | T_c [N m] | fit in frequenza |
|---|---|---|---|
| azimut | 1571 (std 1.9 %) | 0.232 | 75 %, 1.6 dB rms |
| polare | 18 782 (std 7.2 %) | 0.407 | 78 / 77 / 80 % (prima 60 / 69 / 65 %) |

**Anello chiuso** (riferimento di posizione misurato, modello AKD2G, spostamenti da 0.5° a 8°):

- fit di posizione 96–100 % in azimut e 95–99 % in polare, esclusa una finestra anomala;
- picchi di corrente entro il 10–15 % circa;
- si stimano tre guadagni efficaci per asse su un solo doppietto; la corrente non entra nella stima.

Non sono riprodotti l'oscillazione di corrente a ~7 Hz in azimut e lo stick-slip sui piccoli movimenti. In anello
aperto, con la corrente misurata in ingresso, il modello deriva a riposo: la coppia di mantenimento supera T_c.

## Risultati di `run_all.py`

Nella `results_dir` della configurazione:

| percorso | contenuto |
|---|---|
| `report.md` | parametri con i segnaposto evidenziati, validazione della riduzione, frequenze e zeri, simulazioni, stato delle domande aperte |
| `latex/multidof_fem_model.pdf` | documento in inglese con l'analisi FEM in formule e la geometria in TikZ |
| `tables/element_breakdown_<asse>.csv` | scomposizione di ogni elemento per strato (m, J, I_t, G J_p, k) |
| `tables/reduction_validation.csv` | frequenze completo/ridotto/Craig–Bampton e frequenza limite di validità |
| `tables/full_model_frequencies.csv`, `poles_zeros.csv`, `operating_points.csv` | modi, poli, zeri invarianti, equilibri e coppie nominali |
| `tables/modal_modes_<asse>.csv` | tutti i modi del FEM completo: f, ζ, massa modale, c_e, c_m, c_L/c_P, norma di Hankel, tenuto sì/no, controllo di maglia |
| `tables/modal_spectrum.csv` | spettro dei due assi insieme, rigidi per primi |
| `tables/modal_vs_full.csv` | errori in banda contro il modello a tutti i GdL (modale e Guyan), zeri, flessibilità residua, modi aggiunti dalla validazione |
| `tf/ss_th<θ0>_om<Ω0>.{npz,mat,json}` | spazio degli stati del modello principale (in forma modale: anche Φ, ω, ζ e lo stato a blocchi di Gawronski), zeri, num/den |
| `tables/lumped_parameters.csv`, `lumped_vs_full.csv` | modello a parametri concentrati: elementi equivalenti, inerzie nodali, confronto con il modello completo |
| `tables/experimental_comparison.csv` | risonanze e antirisonanze misurate contro il modello concentrato, scarti di livello |
| `tf_concentrato/`, `simulatore_concentrato/` | spazio degli stati del modello concentrato; simulatore autonomo `ne_concentrato.py` / `ne_concentrato.m` |
| `figures/` | Bode 2×2 coppia → encoder, confronti completo/ridotto, forme modali, sensibilità, errore di misura, simulazioni S1–S3, disegno della geometria |

## Come inserire la geometria

Ogni elemento (`id` da 1, nodo `i` → nodo `i+1`) è descritto in uno di due modi:

* `mode: geometry`: `segments` in **serie**, ciascuno con `length_mm` e `layers`
  **coassiali** in parallelo. Ogni strato è `shape: solid` (serve `d_out_mm`) o
  `shape: hollow` (servono `d_in_mm < d_out_mm`), con `material` e
  `carries_torque` (false per distanziali o mozzi che aggiungono solo inerzia).
  Esempio: un elemento fatto da un tratto cavo lungo *a* e uno pieno lungo *b*
  sono due voci di `segments`.
* `mode: manual`: `J_kgm2`, `k_Nm_rad`, `length_mm` dal CAD; per il polare anche
  `m_kg`, `I_t_kgm2`, `y_c_mm` (contributo all'inerzia di azimut).

### Sistemi di riferimento (convenzione De Donatis)

Sezione `frames` di `config.yaml`:

* **N inerziale**: origine O_N (origine dell'assieme Onshape, base dell'albero di
  azimut), versori x̂, ŷ, ẑ (assi del CAD, ẑ verticale lungo l'asse di azimut);
* **corpo 1 (azimut)**: origine o₁ sull'asse di azimut, versori î, ĵ, k̂ = ẑ;
* **corpo 2 (polare)**: origine o₂, centro di rotazione dell'albero polare, versori
  K̂ (puntamento, verso le ottiche), Î (verso l'alto a θ = 0), Ĵ = K̂ × Î = −ĵ (asse
  polare); θ > 0 alza il puntamento.

O_N, o₁ e o₂ stanno su ẑ: `o1_z_mm` = dist(O_N, o₁), `o2_z_mm` = dist(O_N, o₂).
Le origini dei corpi stanno sugli assi di rotazione, non nei baricentri, quindi
l'eccentricità del baricentro è un dato esplicito (d₂ per il telescopio). La posa
CAD è data da `pointing_cad` (K̂ nella posa CAD), `polar_axis_cad` e `theta_cad_deg`.

### Telescopio

Il telescopio **non ha albero**: l'elemento polare 3 è in `mode: telescope` (solo
`length_mm`, per posizioni e disegno). Il telescopio rigido collega i nodi 3 e 4, che
coincidono con P, e la sua inerzia entra solo dalla sezione `payload`:

* `frame: cad_cg` (consigliato): come l'esportazione *Mass properties* di Onshape
  **senza mate connector**, cioè tensore `I_kgm2` rispetto al **baricentro** con assi
  paralleli a x̂, ŷ, ẑ, e baricentro `r_G_mm` misurato da O_N. Il codice lo trasporta
  in o₂ e lo ruota nella terna 2; l'export va fatto nella posa CAD di `frames`.
* `frame: body2_o2`: tensore già rispetto a o₂ nella terna 2, con l'eccentricità
  `d_mm` in terna 2.

`products_sign` fissa la convenzione dei prodotti d'inerzia (`tensor`:
−∫x_i x_j dm); con il solo telescopio `subtract_shaft_parts: false`. La sezione
*Reference frames* del PDF riporta terne, matrici della posa CAD, la corrispondenza
fra componenti in N e nella terna 2 e una checklist per l'esportazione.

### Telescopio sull'azimut (`azimuth.telescope_attachment`, H7)

Con `enabled: true` (serve `polar.hub.enabled: true`) il load di azimut non è più un
corpo unico: forcella + alberi polari + mozzo + anello ruotano con il nodo L (9) e il
telescopio è un nodo P appeso a L dalla molla `k_Nm_rad` attorno a ẑ (catena
e–m–L–p, modello concentrato a 9 GdL, variante `mozzo_azimut`). Solo il modello
concentrato (confronto sperimentale, simulatore, `closedLoopModel`) gestisce il nodo;
l'analisi è `analisi_nodo_azimut.py`.

### Altri campi

Asse polare: `y_node1_from_o2_mm` è la coordinata con segno del nodo 1 lungo Ĵ
(da o₂) e `nodes_increase_along` (`"+J_hat"` o `"-J_hat"`) il verso in cui cresce la
numerazione; nella posa CAD attuale +Ĵ = +x̂. Servono solo per il termine di
Steiner in J^ax.

`k_joint_Nm_rad` aggiunge una cedevolezza di giunzione in serie; `extra_inertias`
aggiunge inerzie concentrate divise a metà sui due nodi. Le inerzie ai nodi
(rotori, forcella) stanno in `node_inertias`: con `J_kgm2`, oppure (solo azimut)
con `frame: cad_cg`, `m_kg`, `r_G_mm`, `I_kgm2` come l'export Onshape, e il codice
usa J = I_zz,G + m d², con d la distanza del baricentro dall'asse di azimut. Se l'export
della forcella comprende gli alberi polari, `includes_polar_shafts: true` evita di
aggiungerli una seconda volta (J^ax non viene sommato).

### Misure sperimentali (`experimental`)

Elenco delle esportazioni CSV di Kollmorgen WorkBench da confrontare con il modello
concentrato (asse, elevazione della prova, telescopio montato o no). La cartella
`Validation` non si usa. `torque_constant_Nm_per_A` converte la coppia del modello nella
corrente della misura. `../LucyIdentification/plot_bode.py --all` traccia le misure da sole,
nello stesso stile delle figure del report.

### Riduzione (`reduction`)

`method: modal` (predefinito) costruisce per ogni asse la base modale del FEM
completo con il telescopio congelato in `theta_ref_deg` e tiene il modo rigido più i
modi flessibili con f ≤ `cutoff_factor` · `f_band_Hz` (ed eventualmente quelli con
norma di Hankel ≥ `hankel_min_ratio` · massima, o un numero forzato `n_flex`). Un
ciclo di validazione confronta il modello lineare troncato con quello a tutti i GdL
in ogni punto di lavoro e aggiunge modi, in ordine di norma di Hankel, finché i
criteri (1 dB, 5°, zeri entro 1 %, 2×2 entro 12 %) sono rispettati.
`static_correction: residual_vectors` aggiunge i vettori statici residui ai nodi
motore e d'aggancio. `method: guyan` e `craig_bampton` restano disponibili. Togliere `placeholder: true` man mano
che i dati diventano reali.

## Codice sorgente

`src/gimbal_mdof/`:

| gruppo | modulo | contenuto |
|---|---|---|
| ingresso | `params.py` | lettura e validazione delle configurazioni, conversione in SI, avvisi |
| FEM | `geometry.py` | proprietà degli elementi da cilindri pieni e cavi |
| | `fem.py` | matrici concentrate e consistenti, nodo P, assemblaggio |
| | `mozzo.py` | mozzo polare e attacco cedevole del telescopio |
| | `payload.py` | telescopio (corpo 2) in o₂ e terna 2; J^ax |
| riduzione | `reduction.py` | Guyan, forma chiusa della catena, Craig–Bampton |
| | `modal.py` | troncamento alla Gawronski, modello su base generica, validazione |
| accoppiamento | `coupled.py` | Lagrange con SymPy sul sistema ridotto |
| | `linearize.py` | equilibrio stazionario e matrici linearizzate |
| | `tf.py` | FRF da P(jω)⁻¹, zeri di Rosenbrock, FdT simboliche, esportazione |
| | `simulate.py` | scenari S1, S2, S3 |
| concentrato | `concentrato.py` | modello concentrato e Newton–Eulero, verificato contro Lagrange |
| | `concentrato_analisi.py` | FdT simboliche e numeriche |
| | `simulatore_concentrato.py` | simulatore non lineare, scenario S4, esportazione Python e MATLAB/Octave |
| misure | `sperimentale.py` | lettura dei CSV WorkBench e confronto col modello, in (rad/s)/A |
| | `vincolo_terra.py` | rigidezza equivalente verso terra k(1+jη) |
| forcella | `modi_forcella.py` | masse modali sullo statore, modello data-driven |
| | `forcella_modale.py` | forma modale, stima su Re/Im con vincoli |
| | `forcella_lugre.py` | LuGre e Dahl con funzione descrittiva, rango 1 modale |
| | `doppietti_tempo.py` | simulatore nel tempo con ingresso corrente, ZOH esatto |
| uscite | `analysis.py` | validazione della riduzione, poli e zeri, sweep |
| | `plotting.py` | figure |
| | `modal_output.py`, `concentrato_output.py`, `sperimentale_output.py` | tabelle e sezioni dei report |
| | `latex_report.py`, `report.py` | documento LaTeX con TikZ, `report.md` |

## Test

```bash
.venv/bin/python -m pytest
```

143 test, tutti superati al 2026-11-01, più 25 test in `closedLoopModel/tests`.

Usano `tests/config_riferimento.yaml`, una configurazione congelata e indipendente da `config.yaml`: modificare i dati
del progetto non li rompe. I `RuntimeWarning` sono trattati come errori. `tests/confronto_guyan_modale.py` è un
prototipo di confronto, non una suite pytest.

| file | test | copertura |
|---|---|---|
| `test_geometry.py` | 12 | proprietà degli elementi e validazione |
| `test_fem.py` | 23 | FEM di ciascun asse |
| `test_reduction.py` | 21 | riduzione di Guyan |
| `test_coupled.py` | 15 | modello accoppiato e linearizzazione |
| `test_payload.py` | 15 | telescopio da configurazione |
| `test_modal.py` | 23 | troncamento modale, accettazione T1–T9 |
| `test_concentrato.py` | 8 | modello concentrato, Newton–Eulero, FdT, simulatore |
| `test_mozzo.py` | 11 | modello hub |
| `test_sperimentale.py` | 4 | lettura dei CSV WorkBench, su dati sintetici |
| `test_vincolo_terra.py` | 3 | vincolo a terra |
| `test_fit.py` | 2 | definizione del fit NRMSE |
| `test_simulate.py` | 2 | scenari S1 e S3 |
| `test_latex.py` | 4 | compilazione del documento LaTeX |

`forcella_modale.py`, `forcella_lugre.py`, `doppietti_tempo.py` e `modi_forcella.py` non hanno test diretti.

## Versionamento

`.gitignore` esclude ambiente virtuale, cache, `results/` e le figure:

- tutti i `*.png` e tutte le cartelle `figures/`;
- i PDF vettoriali di figure salvati fuori da `figures/` (elencati nel file).

I PDF dei report restano versionabili e contengono le figure. Dopo un clone, però, i report non si ricompilano finché
le figure non vengono rigenerate, e rigenerarle richiede `../LucyIdentification`.

## Limiti noti

- **Montaggio.** La risposta sopra i 40 Hz dipende da numero di viti, giochi e precarico: i parametri valgono per una
  configurazione del montaggio.
- **Azimut sotto 10 Hz.** Il comportamento non è ripetibile fra le sessioni e dipende dall'ampiezza.
- **Dahl.** Ha un solo stato e non rappresenta lo stick-slip: la coppia di mantenimento misurata supera T_c.
- **Ampiezza nella funzione descrittiva.** È presa dalla misura, con una circolarità parziale; il bilancio armonico
  autoconsistente non è implementato.
- **Scala del canale AMPS.** Si assume AMPS = 2·I_rms, ma il fattore non è documentato dal PLC.
- **Configurazione dei drive.** I guadagni dei `.akd` non riproducono le sequenze del PLC. I limiti `MOTOR.ICONT` e
  `MOTOR.IPEAK` sembrano impostati con i valori A_dc del datasheet: da verificare sul manuale AKD2G.
- **Stima modale in azimut.** g_Kt = 0.90, inerzia del carico e k₂ finiscono al limite dei vincoli.

## Decisioni e scostamenti dalla specifica

* **Decisioni dell'utente (2026-10-12):** geometria segnaposto; inerzia della
  forcella da input manuale (Q9); materiali del file `.m` (acciaio 8000 kg/m³,
  77 GPa; alluminio 7075 2810 kg/m³, 26.9 GPa); Guyan con Craig–Bampton di verifica.
* **Decisioni dell'utente (2026-10-13):** troncamento modale alla Gawronski come
  riduzione predefinita (Guyan resta come confronto); load della forcella da Onshape
  ("Tavola Rotante [DEF]", alberi polari compresi). Elenco delle modifiche in
  `docs/CAMBIAMENTI_modale.md`.
* **Q10 risolta:** `z_O_cm = 54.9` è già misurato e confermato nel modello a 2 GdL.
* **Ordinamento ω^ℓ ≤ ω^c (sez. 7.4):** vale solo per elementi con inerzia
  simmetrica rispetto alla mezzeria, perché det(M_e^ℓ − M_e^c) = −δ_e², con
  δ_e = (C_e − A_e)/2. Con elementi a più tratti diversi può fallire, e il codice lo
  segnala con un avviso invece di presumerlo. La divisione a metà dell'inerzia
  resta quella richiesta.
* **`bearing_b` dell'azimut** ha chiavi di *elemento*, come il polare e come la
  sezione 7.3 (lo schema della sezione 5 diceva "nodo").
* **Punti di lavoro:** `theta0_deg = [-45, -15, 0, 15, 45]`, cioè il campo
  confermato per il 2 GdL; il default della specifica, `[0, -30, -60, -90]`,
  ne uscirebbe.
* **Test 16 (limite rigido):** con rigidezze ×10⁶ la doppia precisione impone un
  pavimento numerico (FRF confrontate in 1–10 Hz, poli a 10⁻⁴); la spiegazione è nel
  test.
* `gimbal_2dof` ha ricevuto un `pyproject.toml` minimale, necessario per
  l'installazione editable.
