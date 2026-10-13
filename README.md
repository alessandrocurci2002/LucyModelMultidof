# Modello multi-GdL del gimbal alt-azimutale

FEM torsionale degli alberi di azimut e polare, riduzione di Guyan e accoppiamento
dei due assi ridotti tramite l'energia cinetica dell'intero sistema. Implementa
`SPEC_modello_gimbal_multidof.md`; è il secondo modello della catena, dopo il
modello rigido `../gimbal_2dof`, che viene **riusato** (dati CAD del telescopio,
preprocessing inerziale, convenzioni, stile delle figure).

> **Tutti i dati geometrici, i materiali, le inerzie di rotori e forcella e gli
> smorzamenti in `config.yaml` sono segnaposto.** Ogni frequenza, FdT e figura
> prodotta oggi è un *risultato da parametri segnaposto* ed è marcata come tale.

## Installazione

```bash
cd sim/gimbal_multidof
python3 -m venv .venv
.venv/bin/pip install -e ../gimbal_2dof     # modello a 2 GdL, editable
.venv/bin/pip install -e ".[dev]"
```

## Esecuzione

```bash
.venv/bin/python run_all.py                  # tutto: tabelle, figure, simulazioni, LaTeX, report
.venv/bin/python run_all.py --quick          # griglie ridotte, simulazioni brevi (~1 min)
.venv/bin/python run_all.py --no-latex
.venv/bin/python run_all.py --theta0 -15 15 --omega0 0 10 --mass-model consistent --mass-form guyan
.venv/bin/python -m pytest                   # 91 test (sezione 13 della specifica, payload, terne, telescopio)
```

I test usano `tests/config_riferimento.yaml`, una configurazione congelata e
indipendente da `config.yaml`: modificare i dati del progetto non li rompe.

## Risultati (`results/`, rigenerabili)

| percorso | contenuto |
|---|---|
| `report.md` | parametri con i segnaposto evidenziati, validazione della riduzione, frequenze e zeri, simulazioni, stato delle domande aperte |
| `latex/multidof_fem_model.pdf` | documento in inglese con l'analisi FEM in formule e la geometria in TikZ |
| `tables/element_breakdown_<asse>.csv` | scomposizione di ogni elemento per strato (m, J, I_t, G J_p, k) |
| `tables/reduction_validation.csv` | frequenze completo/ridotto/Craig–Bampton e frequenza limite di validità |
| `tables/full_model_frequencies.csv`, `poles_zeros.csv`, `operating_points.csv` | modi, poli, zeri invarianti, equilibri e coppie nominali |
| `tf/ss_th<θ0>_om<Ω0>.{npz,mat,json}` | spazio degli stati (anche per MATLAB), zeri, num/den |
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
  Î (puntamento), Ĵ = ĵ (asse polare), K̂ = Î × Ĵ.

O_N, o₁ e o₂ stanno su ẑ: `o1_z_mm` = dist(O_N, o₁), `o2_z_mm` = dist(O_N, o₂).
Le origini dei corpi stanno sugli assi di rotazione, non nei baricentri, quindi
l'eccentricità del baricentro è un dato esplicito (d₂ per il telescopio). La posa
CAD è data da `pointing_cad` (Î nella posa CAD), `polar_axis_cad` e `theta_cad_deg`.

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

### Altri campi

Asse polare: `y_node1_from_o2_mm` è la coordinata con segno del nodo 1 lungo Ĵ
(da o₂) e `nodes_increase_along` (`"+J_hat"` o `"-J_hat"`) il verso in cui cresce la
numerazione; nella posa CAD attuale +Ĵ = +x̂. Servono solo per il termine di
Steiner in J^ax.

`k_joint_Nm_rad` aggiunge una cedevolezza di giunzione in serie; `extra_inertias`
aggiunge inerzie concentrate divise a metà sui due nodi. Le inerzie ai nodi
(rotori, forcella) stanno in `node_inertias`. Togliere `placeholder: true` man mano
che i dati diventano reali.

## Struttura

```
src/gimbal_mdof/
  params.py        lettura e validazione di config.yaml, conversione SI, avvisi
  geometry.py      proprieta' degli elementi da cilindri pieni e cavi (sez. 6)
  fem.py           matrici concentrate/consistenti, nodo P, assemblaggio (sez. 7)
  reduction.py     Guyan, forma chiusa della catena, Craig-Bampton (sez. 8)
  payload.py       telescopio (corpo 2) da config.yaml, in o_2 e terna 2; J^ax
  coupled.py       Lagrange con SymPy sul sistema ridotto a 7 GdL (sez. 9)
  linearize.py     equilibrio stazionario e matrici linearizzate (sez. 10.1-2)
  tf.py            FRF da P(jw)^-1, zeri di Rosenbrock, FdT simboliche, export
  analysis.py      validazione della riduzione, poli/zeri, sweep
  simulate.py      scenari S1, S2, S3
  plotting.py      figure (stile di gimbal.plotting)
  latex_report.py  documento LaTeX con la geometria in TikZ
  report.py        results/report.md
```

## Decisioni e scostamenti dalla specifica

* **Decisioni dell'utente (2026-10-12):** geometria segnaposto; inerzia della sola
  forcella da input manuale (Q9); materiali del file `.m` (acciaio 8000 kg/m³,
  77 GPa; alluminio 7075 2810 kg/m³, 26.9 GPa); Guyan con Craig–Bampton di verifica.
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
