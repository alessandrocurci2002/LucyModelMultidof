# Specifica per Claude Code — Modello multi-GdL (FEM torsionale + riduzione + accoppiamento) di un telescopio alt-azimutale

## 0. Contesto e obiettivo

Questo è il secondo modello della catena di complessità crescente della tesi. Il primo (`gimbal_2dof/`, specifica `SPEC_modello_gimbal_2dof.md`) tratta forcella e telescopio come due corpi rigidi. Questo modello aggiunge la **flessibilità torsionale degli alberi** di azimut e di elevazione (polare), per ottenere modi strutturali, risonanze e antirisonanze reali.

Il flusso del modello è in quattro stadi, e il codice deve rispecchiarli uno per uno:

1. **Proprietà degli elementi** (inerzia e rigidezza di ogni tratto d'albero), da geometria per composizione di cilindri pieni e cavi oppure da input manuale (CAD).
2. **Modello FEM torsionale** di ciascun asse, isolato: azimut a 9 nodi, polare a 11 nodi, in due varianti di massa (concentrata e consistente).
3. **Riduzione** di ciascun asse a un sistema equivalente a pochi GdL con topologia a catena:
   - azimut: *inerzia encoder – molla – inerzia motore (coppia) – molla – load (forcella)*;
   - polare: *inerzia encoder – molla – inerzia telescopio – molla – inerzia motore (coppia) – molla – inerzia di estremità*.
4. **Accoppiamento** dei due assi ridotti, ricavato dall'energia cinetica dell'intero sistema con la composizione corretta delle rotazioni, esattamente come nel modello a 2 GdL.

**Uscita principale:** la matrice di trasferimento 2×2 dalle coppie $[\tau_\psi,\ \tau_\theta]$ alle posizioni lette dagli encoder $[\psi_e,\ \theta_e]$, nel punto di lavoro $(\theta_0, \Omega_0)$. Inoltre: un documento LaTeX, generato e compilato, che riporta in formule l'analisi FEM svolta e include una rappresentazione grafica della geometria inserita.

**Nota di realismo, da rispettare nei risultati.** Tutti i dati geometrici, i materiali, le rigidezze di giunzione e gli smorzamenti sono oggi **segnaposto**. Finché lo sono, ogni frequenza, FdT o figura che ne dipende deve riportare un avviso visibile ("risultato da parametri segnaposto"), e il report non deve mai presentarla come caratterizzazione del sistema reale.

## 1. Modo di lavorare

- Leggi tutta la specifica prima di scrivere codice, poi proponimi un piano breve.
- **La sezione 3 elenca domande aperte. Prima di implementare i moduli che ne dipendono, chiedimi di risponderle.** Se non rispondo, usa i default indicati, ma ogni output che ne dipende deve riportare un avviso visibile.
- Le convenzioni della sezione 2 sono vincolanti.
- Se trovi incoerenze nei dati, segnalale invece di correggerle in silenzio.
- Non inventare dati (geometrie, materiali, smorzamenti, inerzie di rotori). Se servono e mancano, usa segnaposto dichiarati come tali.
- Implementa un modulo alla volta ed esegui i test della sezione 13 man mano.
- **Riusa** il pacchetto `gimbal_2dof` (dati CAD, preprocessing inerziale del telescopio, convenzioni) importandolo, non duplicandolo: installalo in modalità editable (`pip install -e ../gimbal_2dof`) e dichiaralo in `requirements.txt` / `pyproject.toml`.
- Commenti e docstring in italiano, nomi di variabili coerenti con la notazione qui sotto. Il codice verrà citato in tesi: deve essere leggibile, rigoroso e riproducibile, e ogni risultato tracciabile alle equazioni di questa specifica.

## 2. Convenzioni (vincolanti)

**Ereditate dal modello a 2 GdL, invariate.** Terne W, N, A, B; $\mathbf R_{NB} = \mathbf R_z(\psi)\,\mathbf R_y(\theta)$; θ positivo secondo la mano destra attorno a $y_B$, quindi **θ > 0 abbassa il puntamento** (θ = 0 orizzontale, θ = −90° zenit, ε = −θ solo nelle etichette); ψ positivo antiorario visto dall'alto; tensori in convenzione tensoriale; SI internamente, gradi e Hz in grafici e tabelle.

**Nuove convenzioni di questo modello.**

- **Numerazione.** Nodi numerati da 1 (base 1 nei file di configurazione e nei documenti, base 0 solo internamente nel codice). L'elemento $i$ collega il nodo $i$ al nodo $i+1$. Azimut: nodi 1…9 (elementi 1…8), con il nodo 1 in basso (lato encoder) e il nodo 9 in alto (forcella). Polare: nodi 1…11 (elementi 1…10).
- **Coordinate.** $\psi_i$ è la rotazione assoluta del nodo $i$ dell'albero di azimut attorno a Z. $\theta_i$ è la rotazione del nodo $i$ dell'albero polare **relativa alla forcella**, attorno a $y_A$. Poiché la forcella ruota solo attorno a Z, $\theta_i$ coincide anche con la rotazione assoluta attorno a $y_A$.
- **Encoder.** L'encoder di azimut misura $\psi$ del suo nodo rispetto a terra. L'encoder polare misura $\theta$ del suo nodo rispetto alla forcella (il suo statore è sulla forcella). Entrambi sono quindi direttamente le coordinate generalizzate dei rispettivi nodi.
- **Coppie.** $\tau_\psi$ agisce tra base fissa (statore) e nodo motore di azimut; $\tau_\theta$ tra forcella (statore) e nodo motore polare. Sono positive nel verso delle rispettive coordinate.
- **Unità nei file di input.** Lunghezze e diametri in **mm**; inerzie manuali in **kg·m²**; masse in kg; rigidezze in **N·m/rad**; densità in kg/m³; moduli in Pa. La conversione in SI avviene in un solo punto (`params.py`).

## 3. Decisioni prese, ipotesi e domande aperte

### 3.1 Decisioni confermate dall'utente (vincolanti)

| ID | Decisione | Conseguenza nel modello |
|---|---|---|
| **D1** | **Il telescopio è montato sull'elemento polare 3, tra il nodo 3 e il nodo 4. Si introduce un nodo fittizio P nella mezzeria del telescopio.** | L'elemento 3 viene diviso in P (sezione 7.5). Il telescopio è un corpo rigido collegato a P, e P è un nodo master. L'angolo $\theta_p$ di P è **l'angolo di puntamento del telescopio**: è l'uscita "vera" da confrontare con la lettura dell'encoder $\theta_e$. |
| **D2** | **Gli elementi polari 2 (nodi 2–3) e 4 (nodi 4–5) sono alloggiati sui cuscinetti**, montati nella forcella. | Quegli elementi ruotano in θ come il resto dell'albero e sono portati in azimut dalla forcella. I cuscinetti possono introdurre un attrito viscoso verso la forcella (`bearing_b`, sezione 5). La loro rigidezza radiale non entra nel modello torsionale. Il telescopio si trova quindi tra i due cuscinetti. |
| **D3** | **Il polare ridotto è: encoder con inerzia – molla – inerzia (telescopio) – molla – motore con inerzia – molla – inerzia.** | Nodi master polari: **[1, P, 7, 11]**. GdL ridotti: $\theta_e$ (nodo 1, encoder), $\theta_p$ (nodo P, telescopio), $\theta_m$ (nodo 7, motore), $\theta_r$ (nodo 11, inerzia di estremità). |

### 3.2 Domande aperte e default adottati

Chiedimi di confermare ciascun punto prima di implementare il modulo che ne dipende. Per ogni default, il report deve indicare che è un'ipotesi.

| ID | Domanda | Default adottato se non rispondo | Impatto |
|---|---|---|---|
| **Q1** | **Posizione della mezzeria del telescopio lungo l'elemento 3.** P va nella mezzeria del telescopio (D1). Se il telescopio occupa tutto l'elemento 3, la mezzeria è a metà elemento; altrimenti serve la posizione del suo centro rispetto al nodo 3. | $\xi_p = 0.5$ (telescopio esteso su tutto l'elemento 3). In alternativa, `payload_center_from_node3_mm` → $\xi_p = x_c / L_3$. | Divide la rigidezza dell'elemento 3 tra lato encoder e lato motore. |
| **Q2** | **Il morsetto del telescopio irrigidisce il tratto d'albero che abbraccia?** | No: il tratto conserva la rigidezza dell'albero nudo (ipotesi conservativa). Prevedi `payload_clamp_rigid_mm` (lunghezza del tratto reso rigido attorno a P, default 0). | Rigidezza delle molle adiacenti al telescopio. |
| **Q3** | **Dati del telescopio:** il tensore del corpo B del CAD include l'albero polare (`L01 Albero motore`), che ora è negli elementi. | Tensore e offset del baricentro dal preprocessing del corpo B del modello a 2 GdL, con sottrazione dei pezzi d'albero tramite `shared_component`; avviso finché non è configurata. | Doppio conteggio dell'inerzia assiale. |
| **Q4** | **Metodo di riduzione.** "Riduzione modale" può indicare la condensazione statica di Guyan sui nodi master, Craig–Bampton, o una taratura dei parametri ridotti sulle prime frequenze del modello completo. | Guyan sui nodi master (sezione 8), con Craig–Bampton solo come verifica. | Cambia i valori ridotti; la topologia a catena resta. |
| **Q5** | **Rigidezze:** si calcolano dalla geometria (con G del materiale) o si inseriscono a mano? Le giunzioni (calettamenti, chiavette, giunti, flange) hanno cedevolezze note? | Calcolate dalla geometria dei soli strati portanti, più una cedevolezza di giunzione opzionale in serie (sezione 6). | Le giunzioni dominano spesso la rigidezza reale degli assi a presa diretta. |
| **Q6** | **Materiali** di ogni pezzo (densità ρ e modulo di taglio G). | Acciaio: ρ = 7850 kg/m³, G = 80 GPa. Alluminio: ρ = 2700 kg/m³, G = 26 GPa. Entrambi segnaposto. | Scala lineare su inerzie e rigidezze. |
| **Q7** | **Posizione assiale della catena polare rispetto all'asse di azimut:** distanza con segno del nodo polare 1 dal punto O, lungo $y_A$. | `null` → avviso, e il contributo dell'albero polare all'inerzia di azimut si calcola con i soli momenti trasversali propri (senza Steiner). | Sottostima di $J_\psi$ del load. |
| **Q8** | **Smorzamento.** Non è stato fornito, ma serve per FdT con picchi finiti. | Smorzamento proporzionale alla rigidezza ($\mathbf C = \beta\mathbf K$), con β tarato per ζ = 0.01 sul primo modo flessibile di ciascun asse; attrito dei cuscinetti nullo. Tutto segnaposto. | Ampiezza dei picchi e delle valli delle FRF. |
| **Q9** | **Inerzia del load di azimut (sola forcella).** Il corpo A del CAD include l'albero di azimut, che ora è negli elementi FEM. | Input manuale dell'inerzia al nodo 9 obbligatorio, da un'esportazione Onshape della sola forcella (senza albero di azimut, senza albero polare, senza telescopio). | Doppio conteggio delle inerzie. |
| **Q10** | **`z_O_cm`** (quota del punto O) è ancora da misurare. | Come nel modello a 2 GdL: telescopio bilanciato ($d_z = 0$) con avviso. | Termine gravitazionale. |
| **Q11** | **Velocità massima di inseguimento** $\Omega_{max}$. | 10 °/s, segnaposto. | Entità dell'accoppiamento giroscopico. |
| **Q12** | **Forcella rigida.** Il modello assume che la forcella sia rigida tra il nodo 9 di azimut e i cuscinetti polari (nessuna torsione o flessione della forcella). | Assunzione adottata e dichiarata. | La flessione della forcella è l'unico canale che darebbe lavoro alla reazione di $\tau_\theta$ (sezione 9.4). |

## 4. Topologia del sistema

### 4.1 Modello FEM completo (per asse, disaccoppiato)

```
AZIMUT (dal basso)
                                    τ_ψ
                                     ↓
(1)──e1──(2)──e2──(3)──e3──(4)──e4──(5)──e5──(6)──e6──(7)──e7──(8)──e8──(9)
encoder                                                                load (forcella)

POLARE (da sinistra a destra)
                        TELESCOPIO                               τ_θ
                             │                                    ↓
(1)──e1──(2)──e2──(3)──e3a──(P)──e3b──(4)──e4──(5)──e5──(6)──e6──(7)──e7──(8)──e8──(9)──e9──(10)──e10──(11)
              ▲                            ▲
encoder  cuscinetto                   cuscinetto                                                     estremità
```

La coppia $\tau_\psi$ entra nel nodo 5 di azimut; la coppia $\tau_\theta$ entra nel **nodo 7** polare, dove si trova il rotore del motore (lo statore è sulla forcella). Il nodo P divide l'elemento 3 nella mezzeria del telescopio (sezione 7.5). La numerazione dei nodi 1…11 resta quella dell'utente; P è un'etichetta aggiuntiva, e internamente l'ordine è 1, 2, 3, P, 4, …, 11 (12 GdL).

### 4.2 Sistema ridotto (nodi master)

| Asse | GdL ridotto | Nodo master | Ruolo |
|---|---|---|---|
| Azimut | $\psi_e$ | 1 | encoder |
| Azimut | $\psi_m$ | 5 | motore, coppia $\tau_\psi$ |
| Azimut | $\psi_L$ | 9 | load (forcella), porta l'asse polare |
| Polare | $\theta_e$ | 1 | encoder |
| Polare | $\theta_p$ | P (su elemento 3) | telescopio |
| Polare | $\theta_m$ | 7 | motore, coppia $\tau_\theta$ |
| Polare | $\theta_r$ | 11 | inerzia di estremità |

I nodi master sono parametri (`masters` in `config.yaml`). Encoder, motore, load e telescopio devono essere sempre nodi master: verifica con un `assert`.

Vettore delle coordinate ridotte: $\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ \theta_r]^{\mathsf T}$, 7 GdL.

Catena polare ridotta, con le rigidezze equivalenti in serie:

$$\theta_e \overset{k_{\theta,1}}{\longleftrightarrow} \theta_p \overset{k_{\theta,2}}{\longleftrightarrow} \theta_m \overset{k_{\theta,3}}{\longleftrightarrow} \theta_r$$

$k_{\theta,1}$ raccoglie gli elementi 1, 2 e 3a; $k_{\theta,2}$ gli elementi 3b, 4, 5 e 6; $k_{\theta,3}$ gli elementi 7…10. Il telescopio sta quindi tra encoder e motore: la coppia del motore lo raggiunge attraverso $k_{\theta,2}$, e l'encoder lo legge attraverso $k_{\theta,1}$.

## 5. Input (`config.yaml`)

Schema minimo; la validazione dello schema avviene in `params.py`, con messaggi d'errore chiari.

```yaml
materials:
  acciaio:   {rho: 7850.0, G: 80.0e9, placeholder: true}
  alluminio: {rho: 2700.0, G: 26.0e9, placeholder: true}

azimuth:
  n_nodes: 9
  encoder_node: 1
  torque_node: 5
  load_node: 9
  masters: [1, 5, 9]
  elements:
    - id: 1
      mode: geometry            # geometry | manual
      segments:                 # in SERIE, ordinati dal nodo i al nodo i+1
        - length_mm: 40.0
          layers:               # COASSIALI, in PARALLELO
            - {shape: solid,  d_out_mm: 30.0, material: acciaio, carries_torque: true}
            - {shape: hollow, d_in_mm: 30.0, d_out_mm: 45.0, material: acciaio, carries_torque: false}  # es. distanziale calettato
        - length_mm: 25.0
          layers:
            - {shape: hollow, d_in_mm: 10.0, d_out_mm: 30.0, material: acciaio, carries_torque: true}
      k_joint_Nm_rad: null      # cedevolezza di giunzione in serie (null = assente)
      extra_inertias: []        # es. [{name: "disco encoder", J_kgm2: 1.2e-4, m_kg: 0.05, I_t_kgm2: 6e-5}]
      placeholder: true
    - id: 2
      mode: manual              # inserimento da CAD
      J_kgm2: 3.1e-3
      k_Nm_rad: 2.5e5
      length_mm: 60.0           # serve per la posizione dei nodi e per il disegno
      placeholder: true
    # ... elementi 3..8
  node_inertias:                # inerzie concentrate ai nodi
    - {node: 9, name: "load (forcella)", J_kgm2: null, source: "CAD Onshape, sola forcella", placeholder: true}
  bearing_b: {}                 # {nodo: b in N·m·s/rad}, attrito viscoso verso terra

polar:
  n_nodes: 11
  encoder_node: 1
  torque_node: 7
  end_node: 11
  payload_element: 3            # D1: telescopio tra nodo 3 e nodo 4
  payload_xi: 0.5               # D1/Q1: nodo P nella mezzeria del telescopio (0.5 = telescopio su tutto l'elemento 3)
  payload_center_from_node3_mm: null   # se dato, sostituisce payload_xi: xi_p = x_c / L_3
  payload_clamp_rigid_mm: 0.0   # Q2
  bearing_elements: [2, 4]      # D2: elementi alloggiati sui cuscinetti
  masters: [1, "P", 7, 11]      # D3
  y_node1_from_O_mm: null       # Q7, distanza con segno lungo y_A
  elements: []                  # stesso schema dell'azimut, 10 elementi
  node_inertias: []             # per il polare servono anche m_kg, I_t_kgm2 (posizione = quella del nodo)
  bearing_b: {2: 0.0, 4: 0.0}   # {elemento: b in N·m·s/rad}, attrito viscoso dei cuscinetti verso la FORCELLA (segnaposto)

payload:
  enabled: true
  source: gimbal_2dof           # tensore I_{B,O}|_B e offset d dal preprocessing del 2 GdL
  subtract_shaft_parts: true    # Q3

reduction:
  method: guyan                 # guyan | craig_bampton (solo verifica)
  n_fixed_interface_modes: 0
  mass_form: guyan              # guyan (tridiagonale) | guyan_lumped (diagonale, catena pura)

fem:
  mass_models: [lumped, consistent]
  consistent_shape: linear      # linear | compliance (vedi 7.2)

damping:
  model: stiffness_proportional
  zeta_first_flexible: 0.01     # Q8, segnaposto

operating_points:
  theta0_deg: [0, -30, -60, -90]
  omega0_deg_s: [0, 10]         # Q11
g: 9.80665

frequency:
  f_min_Hz: 0.1
  f_max_Hz: 5000.0              # esteso in automatico a 2× la frequenza propria più alta
  n_points: 4000

latex:
  language: en                  # coerente con gimbal_lagrange_derivation.tex
  compile: true
```

**Regole di validazione.**

- `shape: solid` richiede `d_in_mm` assente o nullo; `shape: hollow` richiede $0 < d_{in} < d_{out}$. Errore altrimenti.
- Ogni segmento di un elemento `geometry` senza `k_Nm_rad` manuale deve avere almeno uno strato con `carries_torque: true`, altrimenti la rigidezza sarebbe nulla: errore esplicito.
- In modalità `manual` per il **polare**, oltre a `J_kgm2` e `k_Nm_rad` servono `m_kg`, `I_t_kgm2` (momento trasversale nel proprio baricentro) e `y_c_mm` (posizione assiale del baricentro dal nodo $i$). Se mancano, avviso: il contributo all'inerzia di azimut di quell'elemento si pone a zero.
- Il numero di elementi deve essere `n_nodes − 1`, con `id` contigui.
- $0 \le$ `payload_xi` $\le 1$. Gli elementi in `bearing_elements` e `payload_element` devono esistere; l'elemento del telescopio non può coincidere con un elemento su cuscinetto.
- Ogni valore con `placeholder: true` propaga l'avviso a tutti gli output che ne dipendono.

## 6. Proprietà degli elementi (`geometry.py`)

Per uno strato anulare di lunghezza $L$, diametri $D_o$ e $D_i$ ($D_i = 0$ per il cilindro pieno), densità ρ e modulo G:

$$m = \rho\,\frac{\pi}{4}\left(D_o^2 - D_i^2\right)L$$

$$J = \frac{\pi\rho L}{32}\left(D_o^4 - D_i^4\right) = \frac{m}{8}\left(D_o^2 + D_i^2\right) \qquad \text{(inerzia assiale)}$$

$$I_t = m\left[\frac{D_o^2 + D_i^2}{16} + \frac{L^2}{12}\right] \qquad \text{(inerzia trasversale nel proprio baricentro)}$$

$$J_p = \frac{\pi}{32}\left(D_o^4 - D_i^4\right) \qquad \text{(momento polare d'area)}$$

**Composizione.**

- Strati coassiali di un segmento: inerzie e masse si sommano; le rigidezze degli strati portanti si sommano in parallelo, $k_s = \dfrac{1}{L_s}\displaystyle\sum_{\ell\ \text{portanti}} G_\ell J_{p,\ell}$.
- Segmenti di un elemento: inerzie e masse si sommano; le cedevolezze si sommano in serie, $\dfrac{1}{k_e} = \displaystyle\sum_s \frac{1}{k_s} + \frac{1}{k_{joint}}$.
- Inerzie extra (`extra_inertias`): si sommano a $J_e$ e, per il polare, a $m$ e $I_t$.
- Ogni strato ha posizione assiale nota (inizio del segmento + $L_s/2$): serve per la massa consistente (sezione 7.2) e per il contributo all'inerzia di azimut (sezione 9.2).

Ogni funzione restituisce anche una scomposizione per strato (tabella CSV `results/tables/element_breakdown_<asse>.csv`), per rendere tracciabile ogni numero.

## 7. Modello FEM per asse (`fem.py`)

### 7.1 Massa concentrata

L'inerzia di ogni elemento si divide a metà sui suoi due nodi:

$$\mathbf M_e^{\ell} = \frac{J_e}{2}\begin{bmatrix}1&0\\0&1\end{bmatrix}$$

### 7.2 Massa consistente

Con funzioni di forma lineari nella coordinata fisica $x \in [0, L_e]$, $N_1 = 1 - x/L_e$, $N_2 = x/L_e$, e inerzia per unità di lunghezza $\bar J(x)$ **costante a tratti** (un valore per segmento):

$$(\mathbf M_e^{c})_{ab} = \sum_{s}\bar J_s\int_{x_s}^{x_s + L_s} N_a(x)\,N_b(x)\,dx$$

Gli integrali vanno calcolati in forma chiusa, segmento per segmento. Con un solo segmento uniforme si ritrova $\mathbf M_e^{c} = \frac{J_e}{6}\begin{bmatrix}2&1\\1&2\end{bmatrix}$ (test). Le inerzie extra di un elemento (non distribuite) si dividono a metà sui due nodi in entrambi i modelli.

Opzione `consistent_shape: compliance`: invece di funzioni di forma lineari in $x$, usa la forma statica esatta per rigidezza variabile, lineare nella cedevolezza accumulata. Implementala come estensione facoltativa, disattivata di default.

### 7.3 Rigidezza e assemblaggio

$$\mathbf K_e = k_e\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$$

Assembla $\mathbf M$ e $\mathbf K$ globali per somma sui nodi condivisi. Aggiungi le inerzie concentrate (`node_inertias`) sulla diagonale. Il vettore delle forze è $\mathbf Q = \tau\,\mathbf e_p$, con $p$ il nodo motore.

**Telescopio e load nei modelli del singolo asse.** Le analisi di un asse isolato (FEM completo e catena ridotta isolata) servono a validare la riduzione e devono contenere le inerzie dominanti, altrimenti non rappresentano nulla di fisico:

- polare: al nodo P si somma $I_{yy,p}$ (inerzia del telescopio attorno all'asse polare, per O), e alla rigidezza di P si somma $k_g(\theta_0)$;
- azimut: al nodo 9 si somma $J^{ax} + J_{\psi,p}(\theta_0)$ (sezione 9.2).

Nel modello accoppiato, invece, telescopio e $J^{ax}$ entrano solo attraverso $T$ e $V$ (sezione 9) e **non vanno sommati due volte**. Tieni le due costruzioni separate nel codice (un flag `isolated=True/False`).

**Smorzamento.** $\mathbf C = \beta\mathbf K + \mathbf C_b$, con β tale che ζ = `zeta_first_flexible` sul primo modo flessibile ($\beta = 2\zeta/\omega_1$). La parte proporzionale alla massa è esclusa di proposito: smorzerebbe il moto rigido rispetto allo spazio inerziale, cosa che un albero libero non fa.

$\mathbf C_b$ è l'attrito viscoso dei cuscinetti (D2). Per un elemento su cuscinetto con coefficiente $b$, il contributo è $\frac{b}{2}\operatorname{diag}(1, 1)$ sui suoi due nodi. Per il polare il cuscinetto è nella forcella: l'attrito agisce sulla rotazione relativa alla forcella, che è proprio θ. Per l'azimut agisce verso terra.

### 7.4 Controlli

- $\mathbf M$, $\mathbf K$ simmetriche; $\mathbf M$ definita positiva.
- $\mathbf K\mathbf 1 = \mathbf 0$: un solo modo rigido per asse.
- $\mathbf 1^{\mathsf T}\mathbf M\mathbf 1 = \sum_e J_e + \sum J_{nodi}$, in entrambi i modelli di massa.
- $\mathbf M^{\ell} - \mathbf M^{c}$ semidefinita positiva (con funzioni di forma lineari), quindi $\omega_k^{\ell} \le \omega_k^{c}$ per ogni $k$ (dimostrazione nel documento LaTeX).

**Uscite del modulo:** frequenze e forme modali del modello completo di ciascun asse, isolato, per entrambi i modelli di massa; FdT coppia → encoder e coppia → nodo motore (collocata) del modello completo, che servono a validare la riduzione.

### 7.5 Nodo ausiliario P sull'elemento 3 polare (telescopio)

Il telescopio è montato sull'elemento 3 (D1). Perché sia un unico GdL del sistema ridotto, l'elemento 3 viene diviso nel punto di attacco $x_P = \xi_p L_3$ in due sotto-elementi, 3a (nodo 3 → P) e 3b (P → nodo 4):

- **modalità `geometry`:** i segmenti dell'elemento 3 vengono tagliati in $x_P$ (un segmento attraversato da $x_P$ diventa due segmenti con gli stessi strati). Inerzie e rigidezze di 3a e 3b si ricalcolano con la sezione 6. Se `payload_clamp_rigid_mm` > 0, il tratto di quella lunghezza centrato in P ha cedevolezza nulla;
- **modalità `manual`:** si assume l'elemento uniforme, quindi $k_{3a} = k_3/\xi_p$, $k_{3b} = k_3/(1 - \xi_p)$, $J_{3a} = \xi_p J_3$, $J_{3b} = (1 - \xi_p)J_3$, con un avviso.

Con $\xi_p = 0$ o $\xi_p = 1$ non si divide nulla: P coincide con il nodo 3 o con il nodo 4.

Controlli: $1/k_{3a} + 1/k_{3b} = 1/k_3$ e $J_{3a} + J_{3b} = J_3$ (con tratto rigido nullo).

Nelle figure, nelle tabelle e nel LaTeX il nodo si chiama **P** e la numerazione 1…11 dell'utente resta invariata.

## 8. Riduzione (`reduction.py`)

### 8.1 Condensazione di Guyan sui nodi master

Partiziona i GdL in master (m) e slave (s):

$$\mathbf T = \begin{bmatrix}\mathbb 1\\ -\mathbf K_{ss}^{-1}\mathbf K_{sm}\end{bmatrix}, \qquad \mathbf M_r = \mathbf T^{\mathsf T}\mathbf M\,\mathbf T, \quad \mathbf K_r = \mathbf T^{\mathsf T}\mathbf K\,\mathbf T, \quad \mathbf C_r = \mathbf T^{\mathsf T}\mathbf C\,\mathbf T$$

$\mathbf K_{ss}$ è sempre invertibile, purché ci sia almeno un master.

**Forma chiusa per una catena (target di regressione).** Tra due master consecutivi $a < b$, la rigidezza equivalente è quella delle molle in serie, e lo spostamento statico di uno slave $j$ è un'interpolazione pesata sulla cedevolezza:

$$\frac{1}{k_{ab}} = \sum_{i=a}^{b-1}\frac{1}{k_i}, \qquad \varphi_j = (1 - w_j)\,\varphi_a + w_j\,\varphi_b, \qquad w_j = \frac{\sum_{i=a}^{j-1} 1/k_i}{\sum_{i=a}^{b-1} 1/k_i}$$

Gli slave oltre l'ultimo master (o prima del primo) seguono rigidamente il master adiacente. Ne segue che $\mathbf K_r$ è tridiagonale (catena di molle equivalenti) e anche $\mathbf M_r$ è tridiagonale.

### 8.2 Forma della massa ridotta

- `mass_form: guyan`: $\mathbf M_r$ tridiagonale, come esce dalla condensazione (più accurata).
- `mass_form: guyan_lumped`: somma per righe di $\mathbf M_r$ sulla diagonale. Si ottiene la catena pura "inerzia – molla – inerzia" descritta in sezione 0. Conserva l'inerzia totale e il modo rigido.

Calcola sempre entrambe e riporta il confronto.

### 8.3 Validazione della riduzione (obbligatoria, nel report)

- Tabella delle frequenze proprie: modello completo contro ridotto, per ogni combinazione (concentrata/consistente) × (guyan/guyan_lumped), con errore relativo.
- Guyan è un'approssimazione di Rayleigh–Ritz del modello completo, quindi con `mass_form: guyan` deve valere $\omega_k^{r} \ge \omega_k^{full}$ per ogni $k$ (test).
- Sovrapposizione delle FRF coppia → encoder: completo contro ridotto, con indicazione della frequenza oltre la quale il ridotto non è più rappresentativo.
- Craig–Bampton con $n$ modi a interfaccia fissa (`n_fixed_interface_modes`) solo come riferimento di convergenza. Con $n = 0$ deve coincidere con Guyan (test).

## 9. Modello accoppiato (`coupled.py`)

### 9.1 Cinematica

La forcella ruota con $\psi_L$ (ipotesi Q12). Ogni corpo $j$ dell'asse polare ha velocità angolare $\boldsymbol\omega_j = \dot\psi_L\,\mathbf Z + \dot\theta_j\,\mathbf y_A$. Per il payload, in terna B, come nel modello a 2 GdL:

$$\boldsymbol\omega_p\big|_B = \begin{bmatrix}-\dot\psi_L\sin\theta_p\\ \dot\theta_p\\ \dot\psi_L\cos\theta_p\end{bmatrix}$$

Ricavala con SymPy dalla formula di Poisson e verifica che coincida.

### 9.2 Energia cinetica

$$T = \tfrac12\,\dot{\mathbf q}_\psi^{\mathsf T}\mathbf M_{r,\psi}\,\dot{\mathbf q}_\psi \;+\; \tfrac12\,J^{ax}\,\dot\psi_L^2 \;+\; \tfrac12\,\dot{\mathbf q}_\theta^{\mathsf T}\mathbf M_{r,\theta}\,\dot{\mathbf q}_\theta \;+\; \tfrac12\,\boldsymbol\omega_p|_B^{\mathsf T}\,\mathbf I_{p,O}|_B\,\boldsymbol\omega_p|_B$$

con $\mathbf q_\psi = [\psi_e, \psi_m, \psi_L]$ e $\mathbf q_\theta = [\theta_e, \theta_p, \theta_m, \theta_r]$. Il telescopio ruota con $\theta_p$, il GdL del nodo P.

**Termine $J^{ax}$, da documentare in tesi.** Ogni parte assialsimmetrica dell'albero polare ha tensore $\operatorname{diag}(I_t, J, I_t)$ in qualunque terna ruotata attorno a $y_A$. Quindi la sua energia cinetica è

$$\tfrac12\left(I_t + m\,y_c^2\right)\dot\psi_L^2 + \tfrac12\,J\,\dot\theta_j^2$$

senza termini misti e senza dipendenza da θ. Ne seguono due fatti:

1. la parte assiale ($J$) è già in $\mathbf M_{r,\theta}$;
2. la parte trasversale si somma in una costante, $J^{ax} = \sum_{\text{strati e inerzie polari}}\left(I_t + m\,y_c^2\right)$, con $y_c$ la distanza del baricentro dello strato da O lungo $y_A$ (Q7).

**Gli alberi polari non producono quindi accoppiamento giroscopico.** L'accoppiamento dinamico nasce solo dal payload (tensore anisotropo, baricentro fuori asse). Questo è il motivo per cui il payload deve stare su un nodo master. Il report deve riportare questa osservazione.

L'inerzia del load (forcella) è già in $\mathbf M_{r,\psi}$, perché il nodo 9 è master e l'inerzia di un master passa invariata attraverso Guyan. $\mathbf I_{p,O}|_B$ viene dal preprocessing del 2 GdL, depurato dei pezzi d'albero già negli elementi (Q9).

### 9.3 Energia potenziale e dissipazione

$$V = \tfrac12\,\mathbf q_\psi^{\mathsf T}\mathbf K_{r,\psi}\,\mathbf q_\psi + \tfrac12\,\mathbf q_\theta^{\mathsf T}\mathbf K_{r,\theta}\,\mathbf q_\theta + m_p\,g\left(-d_x\sin\theta_p + d_z\cos\theta_p\right)$$

$$\mathcal D = \tfrac12\,\dot{\mathbf q}_\psi^{\mathsf T}\mathbf C_{r,\psi}\,\dot{\mathbf q}_\psi + \tfrac12\,\dot{\mathbf q}_\theta^{\mathsf T}\mathbf C_{r,\theta}\,\dot{\mathbf q}_\theta$$

Gli alberi sono sull'asse: non hanno energia gravitazionale.

### 9.4 Forze generalizzate e coppia di reazione

$$\mathbf Q = \tau_\psi\,\mathbf e_{\psi_m} + \tau_\theta\,\mathbf e_{\theta_m}$$

**Derivazione da riportare nel codice e nel documento.** Il motore polare applica $+\tau_\theta\,\mathbf y_A$ al rotore (nodo 7) e $-\tau_\theta\,\mathbf y_A$ allo statore (forcella). Le rotazioni virtuali sono $\delta\boldsymbol\phi_{rotore} = \delta\psi_L\,\mathbf Z + \delta\theta_m\,\mathbf y_A$ e $\delta\boldsymbol\phi_{statore} = \delta\psi_L\,\mathbf Z$. Il lavoro virtuale è quindi

$$\delta W = \tau_\theta\,\mathbf y_A\cdot\left(\delta\boldsymbol\phi_{rotore} - \delta\boldsymbol\phi_{statore}\right) = \tau_\theta\,\delta\theta_m$$

La coppia di reazione sul nodo 9 di azimut esiste fisicamente, ma è diretta lungo $\mathbf y_A \perp \mathbf Z$. Sull'unico GdL della forcella ($\psi_L$, rotazione attorno a Z) compie lavoro nullo: la scaricano i cuscinetti di azimut come momento flettente. Con le coordinate polari relative alla forcella, lo scambio di reazione è incluso automaticamente e **non va aggiunto a mano come forza su $\psi_L$**. Comparirebbe nella dinamica solo con un GdL di flessione della forcella, che è fuori scope (Q12).

### 9.5 Equazioni di Lagrange

$$\frac{d}{dt}\frac{\partial\mathcal L}{\partial\dot{\mathbf q}} - \frac{\partial\mathcal L}{\partial\mathbf q} + \frac{\partial\mathcal D}{\partial\dot{\mathbf q}} = \mathbf Q, \qquad \mathcal L = T - V$$

Derivale con SymPy, senza codificarle a mano. Costruisci $\mathbf C(\mathbf q, \dot{\mathbf q})$ con i simboli di Christoffel e verifica che $\dot{\mathbf M} - 2\mathbf C$ sia antisimmetrica. Esporta le equazioni in LaTeX.

**Target di regressione** con `mass_form: guyan_lumped` e prodotti d'inerzia del payload nulli:

$$\left[J_{L}^{*} + J_{\psi,p}(\theta_p)\right]\ddot\psi_L + J_{\psi,p}'(\theta_p)\,\dot\theta_p\,\dot\psi_L - k_{\psi,2}(\psi_m - \psi_L) + \ldots = 0$$

$$\left[M_{r,\theta,pp} + I_{yy,p}\right]\ddot\theta_p - \tfrac12 J_{\psi,p}'(\theta_p)\,\dot\psi_L^2 + \frac{\partial V_g}{\partial\theta_p} - k_{\theta,1}(\theta_e - \theta_p) + k_{\theta,2}(\theta_p - \theta_m) + \ldots = 0$$

$$J_{\psi,p}(\theta) = I_{xx,p}\sin^2\theta + I_{zz,p}\cos^2\theta, \qquad J_L^{*} = M_{r,\psi,LL} + J^{ax}$$

dove $k_{\psi,2}$ è la rigidezza equivalente della molla ridotta che arriva a $\psi_L$, e $k_{\theta,1}$, $k_{\theta,2}$ quelle delle due molle ridotte ai lati del telescopio (sezione 4.2). Le equazioni di $\psi_e, \psi_m, \theta_e, \theta_m, \theta_r$ devono restare lineari e identiche a quelle delle catene ridotte isolate.

## 10. Linearizzazione e funzioni di trasferimento (`linearize.py`, `tf.py`)

### 10.1 Punto di lavoro

Assegnati $\theta_0$ e $\Omega_0$: $\theta_p = \theta_0$, tutte le $\dot\psi = \Omega_0$, tutte le $\dot\theta = 0$, $\ddot{\mathbf q} = 0$. Le altre coordinate (torsioni statiche delle catene dovute alla coppia di mantenimento e all'attrito) e le coppie nominali $\boldsymbol\tau_0$ si ottengono risolvendo l'equilibrio stazionario (sistema lineare più termini del payload, con `scipy.optimize.fsolve` o in forma chiusa). Fissa $\psi_e = 0$ (ψ ciclica). Riporta sempre $\tau_{\theta 0}$ e le torsioni statiche della catena polare.

**Risultato atteso, da verificare e commentare in tesi.** La coppia di mantenimento passa dal motore al telescopio solo attraverso $k_{\theta,2}$. Il ramo encoder ($k_{\theta,1}$) e il ramo di estremità ($k_{\theta,3}$) non trasmettono coppia statica (in assenza di attrito statico). Quindi, in statica, $\theta_e = \theta_p$ e $\theta_r = \theta_m$, mentre $\theta_m - \theta_p = \tau_{\theta 0}/k_{\theta,2}$. L'encoder polare, montato dalla parte del telescopio opposta al motore, legge esattamente l'angolo del telescopio in condizioni stazionarie e se ne discosta solo in dinamica. Lo stesso vale per l'azimut: senza attrito ai cuscinetti tra motore e load, $\psi_e = \psi_m = \psi_L$ a velocità costante.

### 10.2 Matrici linearizzate

Ricavale come jacobiani simbolici del residuo, valutati nel punto di lavoro, esattamente come nel 2 GdL:

$$\mathbf M_0\,\delta\ddot{\mathbf q} + \mathbf C_0\,\delta\dot{\mathbf q} + \mathbf K_0\,\delta\mathbf q = \delta\boldsymbol\tau$$

**Target con prodotti nulli.** Rispetto alle catene ridotte isolate cambiano solo questi termini:

- $\mathbf C_0[\psi_L, \theta_p] = +\gamma$ e $\mathbf C_0[\theta_p, \psi_L] = -\gamma$, con $\gamma = (I_{xx,p} - I_{zz,p})\,\Omega_0\sin 2\theta_0$;
- $\mathbf K_0[\theta_p, \theta_p]$ cresce di $k_c + k_g$, con $k_c = -(I_{xx,p} - I_{zz,p})\,\Omega_0^2\cos 2\theta_0$ e $k_g = m_p\,g\,(d_x\sin\theta_0 - d_z\cos\theta_0)$.

### 10.3 Spazio degli stati e uscite

$\mathbf x = [\delta\mathbf q,\ \delta\dot{\mathbf q}]$ (14 stati), con $\mathbf A$ e $\mathbf B_u$ come nel 2 GdL. Uscite:

| Uscita | Significato | Ruolo |
|---|---|---|
| $\psi_e$, $\theta_e$ | encoder | **uscite principali** |
| $\psi_m$, $\theta_m$ | nodi motore | FdT collocate, di riferimento |
| $\psi_L$, $\theta_p$ | load e telescopio (nodo P) | **puntamento reale**: confrontato con gli encoder, dà l'errore tra angolo misurato e angolo vero |
| $\theta_r$ | inerzia di estremità polare | diagnostica |

**Consegna principale: la matrice 2×2 coppia → encoder.**

$$\begin{bmatrix}\Psi_e\\ \Theta_e\end{bmatrix} = \begin{bmatrix}G_{\psi_e\tau_\psi} & G_{\psi_e\tau_\theta}\\ G_{\theta_e\tau_\psi} & G_{\theta_e\tau_\theta}\end{bmatrix}\begin{bmatrix}T_\psi\\ T_\theta\end{bmatrix}$$

per ogni punto di lavoro, per ogni modello di massa e per ogni forma della massa ridotta.

**Consegna secondaria: errore di misura.** Calcola anche le FdT coppia → puntamento reale ($\psi_L$, $\theta_p$) e coppia → errore di misura ($\psi_e - \psi_L$, $\theta_e - \theta_p$), in un Bode dedicato. Quantificano quanto l'encoder si discosta dall'angolo vero del carico e del telescopio al variare della frequenza; in statica l'errore è nullo (sezione 10.1).

- **FRF:** calcola $\mathbf G(j\omega)$ da $\mathbf P(j\omega)^{-1} = (\mathbf M_0(j\omega)^2 + \mathbf C_0 j\omega + \mathbf K_0)^{-1}$, non dai polinomi.
- **Poli:** autovalori di $\mathbf A$, classificati come nel 2 GdL (nulli, reali, complessi coniugati con $f_n$ e ζ).
- **Zeri di ogni FdT SISO:** zeri invarianti dal pencil di Rosenbrock $\begin{bmatrix}\mathbf A - s\mathbb 1 & \mathbf b_j\\ \mathbf c_i & 0\end{bmatrix}$, con `scipy.linalg.eig` generalizzato. Scarta gli autovalori infiniti.
- **FdT simboliche:** solo per le catene ridotte isolate (3 e 4 GdL, con $\Omega_0 = 0$), esportate in LaTeX. Per il sistema accoppiato bastano le matrici $\mathbf M_0, \mathbf C_0, \mathbf K_0$ simboliche e le FdT numeriche.
- **Esportazione:** `results/tf/` con le matrici di stato in `.npz` e `.mat` (`scipy.io.savemat`, per l'uso in MATLAB) e i coefficienti numeratore/denominatore in JSON.

**Nota da riportare.** La FdT coppia → encoder è **non collocata**: encoder e motore stanno su nodi diversi. L'alternanza poli-zeri non è garantita e possono comparire zeri a parte reale positiva o coppie di zeri mancanti. La FdT coppia → nodo motore è invece collocata e deve mostrare l'alternanza (test).

## 11. Analisi, simulazioni e figure (`analysis.py`, `simulate.py`, `plotting.py`)

1. Tabelle delle frequenze proprie e degli zeri per tutte le combinazioni di modello e punto di lavoro (CSV).
2. Forme modali di ciascun asse: modello completo, con i master evidenziati, e modello ridotto sovrapposto.
3. Bode 2×2 coppia → encoder con poli e zeri segnati; sovrapposizione concentrata/consistente e completo/ridotto ($\Omega_0 = 0$, accoppiamento spento).
4. Sensibilità: frequenze in funzione di $\theta_0$ (una curva per $\Omega_0$); FRF incrociate in funzione di $\Omega_0$.
5. Simulazioni non lineari del sistema accoppiato ridotto (`solve_ivp`, `DOP853`, `rtol=1e-10`, `atol=1e-12`):
   - **S1:** doppietto su $\tau_\psi$ a $\Omega_0 = 0$; confronto con il lineare, errore relativo sotto l'1%.
   - **S2:** inseguimento a $\Omega_{max}$ con doppietto su $\tau_\psi$; risposta incrociata su $\theta_e$ e $\theta_p$.
   - **S3:** conservazione, con $\mathbf C = 0$ e $\boldsymbol\tau = 0$: deriva relativa di $E = T + V$ e di $p_\psi = \sum\partial\mathcal L/\partial\dot\psi_i$ sotto $10^{-8}$.

Figure in PDF vettoriale e PNG a 300 dpi, etichette in italiano con unità; ogni figura da segnaposto ha l'avviso nel titolo.

## 12. Documento LaTeX (`latex_report.py`)

Genera `results/latex/multidof_fem_model.tex` e compilalo con `latexmk -pdf` (oppure `pdflatex` tre volte). Se LaTeX non è installato, lascia il `.tex`, stampa un avviso e salta il test di compilazione. In lingua inglese (`latex.language`), coerente con `gimbal_lagrange_derivation.tex`, che puoi riusare come base per le derivazioni generali.

**Contenuti, nell'ordine:**

1. Descrizione del sistema e convenzioni.
2. **Rappresentazione grafica della geometria inserita**, generata da `config.yaml` in TikZ, una figura per asse. Vista laterale stilizzata in cui ogni segmento è un rettangolo con larghezza proporzionale a $L$ e altezza proporzionale a $D_o$; foro dei cilindri cavi con un riempimento diverso; strati non portanti (distanziali) con un colore distinto. Quote $L$, $D_o$, $D_i$ su ogni segmento, numerazione di nodi ed elementi, frecce delle coppie, encoder e load; per il polare, i cuscinetti sugli elementi 2 e 4 (simbolo di cuscinetto sopra e sotto l'albero) e il telescopio sull'elemento 3, con il nodo P nel punto $\xi_p$. Se le proporzioni reali rendono illeggibili i pezzi piccoli, usa una scala per elemento e dichiaralo in didascalia. Salva anche un PNG del disegno per un controllo rapido.
3. Formule delle proprietà degli elementi (sezione 6), con una tabella per elemento: contributi per strato, $J_e$, $k_e$, $m$, $I_t$.
4. Formulazione FEM: matrici di elemento concentrata e consistente (integrale a tratti della 7.2), assemblaggio, matrici globali simboliche e numeriche, dimostrazione $\mathbf M^\ell \succeq \mathbf M^c \Rightarrow \omega^\ell \le \omega^c$.
5. Riduzione di Guyan: formule, forma chiusa per la catena, parametri ridotti (tabella), confronto delle frequenze.
6. Modello accoppiato: cinematica, $T$ con il termine $J^{ax}$ e la sua giustificazione, $V$, derivazione del lavoro virtuale della coppia di reazione (9.4), equazioni di moto.
7. Linearizzazione, matrici $\mathbf M_0, \mathbf C_0, \mathbf K_0$, FdT simboliche delle catene isolate, Bode coppia → encoder (figure incluse), tabelle di poli e zeri.
8. Elenco delle ipotesi e dei segnaposto (sezione 3), con il loro stato.

**Note tecniche obbligatorie** (problemi già incontrati):

- Nel preambolo `\setcounter{MaxMatrixCols}{20}`: una `bmatrix` di amsmath accetta di default al massimo 10 colonne, e la matrice polare è 11×11.
- Vettori di forzante come $\tau\,\mathbf e_p$, non come vettori riga espliciti (vanno a capo).
- Blocchi di equazioni in `gather` (centrati), non `align` ancorato all'uguale.
- `\texorpdfstring` per la matematica nei titoli di sezione.
- La compilazione deve chiudersi senza errori né warning di riferimenti non risolti; controlla il log in un test.

## 13. Test (`pytest`)

**Geometria**

1. Cilindro cavo con $D_i \to 0$ uguale al pieno; additività (cavo di lunghezza $a$ più pieno di lunghezza $b$ uguale alla somma dei singoli).
2. Valori analitici di $J$, $m$, $I_t$, $k$ per un cilindro noto; due segmenti uguali in serie danno $k/2$; due strati portanti in parallelo danno la somma.
3. Errori di validazione: pieno con $D_i > 0$, cavo con $D_i \ge D_o$, segmento senza strati portanti.

**FEM**

4. Matrici di elemento: con un segmento uniforme, $\mathbf M^c = \frac{J}{6}\begin{bmatrix}2&1\\1&2\end{bmatrix}$ e $\mathbf M^\ell = \frac{J}{2}\mathbb 1$.
5. Simmetria, $\mathbf K\mathbf 1 = \mathbf 0$, conservazione dell'inerzia totale, un solo autovalore nullo.
6. $\omega_k^\ell \le \omega_k^c$ per ogni modo; identità $v^{\mathsf T}(\mathbf M^\ell - \mathbf M^c)v = \sum_i \frac{J_i}{6}(v_i - v_{i+1})^2$ per elementi uniformi.
7. Nodo P: $1/k_{3a} + 1/k_{3b} = 1/k_3$ e $J_{3a} + J_{3b} = J_3$; con $\xi_p \in \{0, 1\}$ nessuna divisione; in modalità `geometry`, un segmento tagliato in due conserva massa e inerzia; con $\xi_p = 0.5$ e massa concentrata, l'inerzia di P vale $J_3/2$ (più $I_{yy,p}$ nel modello isolato).
8. Cuscinetti: con `bearing_b` non nullo su un elemento, $\mathbf C_b$ ha $b/2$ sui suoi due nodi; nell'azimut isolato, uno dei due poli nulli del moto rigido si sposta in $-\sum b / \sum J$.
9. Convergenza: raffinando un albero uniforme, le frequenze tendono a quelle analitiche della barra libera-libera, $\omega_n = \frac{n\pi}{L}\sqrt{G/\rho}$, dal basso con massa concentrata e dall'alto con consistente.

**Riduzione**

10. $\mathbf K_r$ e il vettore dei pesi $w_j$ coincidono con la forma chiusa della 8.1.
11. Cedevolezza statica tra i master identica tra completo e ridotto; inerzia totale conservata in entrambe le forme della massa.
12. Con `mass_form: guyan`: $\omega_k^{r} \ge \omega_k^{full}$. Craig–Bampton con $n = 0$ uguale a Guyan.

**Modello accoppiato**

13. La velocità angolare da Poisson coincide con la 9.1.
14. Con prodotti nulli, le equazioni simboliche coincidono con il target della 9.5; $\dot{\mathbf M} - 2\mathbf C$ antisimmetrica.
15. **Disaccoppiamento:** con $\Omega_0 = 0$ e prodotti del payload nulli, le FdT incrociate sono nulle (entro la tolleranza numerica), e quelle dirette coincidono con le catene ridotte isolate, con $J_L^{*} + J_{\psi,p}(\theta_0)$ su $\psi_L$, $I_{yy,p}$ e $k_g$ su $\theta_p$.
16. **Limite rigido:** moltiplicando tutte le rigidezze per $10^6$, i poli a bassa frequenza e le FRF sotto la prima risonanza coincidono con il modello a 2 GdL (importato da `gimbal_2dof`) con inerzie totali equivalenti.
17. **Collocazione:** con smorzamento nullo, $\Omega_0 = 0$ e asse isolato, la FdT coppia → nodo motore mostra poli e zeri alternati.
18. Coerenza tra le rappresentazioni: FRF da spazio degli stati uguale a FRF da $\mathbf P(j\omega)^{-1}$; poli uguali alle radici di $\det\mathbf P(s)$.
19. Il lavoro virtuale di $\tau_\theta$ non compare nell'equazione di $\psi_L$ (9.4).
20. Scenari S1 e S3 della sezione 11 entro le tolleranze.
21. **Equilibrio statico polare** (sezione 10.1): $\theta_e = \theta_p$, $\theta_r = \theta_m$ e $\theta_m - \theta_p = \tau_{\theta 0}/k_{\theta,2}$, entro la tolleranza numerica.

**Documento**

22. Il `.tex` compila senza errori (saltato se LaTeX non c'è) e contiene una figura di geometria per asse.

## 14. Struttura del progetto e requisiti tecnici

La cartella va creata nella radice del repository, accanto a `gimbal_2dof/`.

```
gimbal_multidof/
├── config.yaml
├── pyproject.toml / requirements.txt   # numpy, scipy, sympy, matplotlib, pyyaml, pytest; gimbal_2dof (editable)
├── README.md                # come eseguire, dove sono i risultati, ipotesi e segnaposto
├── run_all.py               # CLI: --theta0, --omega0, --mass-model, --mass-form, --quick, --no-latex
├── src/gimbal_mdof/
│   ├── params.py            # lettura e validazione di config.yaml, conversione SI, propagazione dei segnaposto
│   ├── geometry.py          # sezione 6
│   ├── fem.py               # sezione 7
│   ├── reduction.py         # sezione 8
│   ├── payload.py           # adattatore verso gimbal_2dof (tensore del telescopio, d, sottrazioni)
│   ├── coupled.py           # sezione 9
│   ├── linearize.py         # sezione 10.1–10.2
│   ├── tf.py                # sezione 10.3: FRF, poli, zeri, esportazione
│   ├── analysis.py
│   ├── simulate.py
│   ├── plotting.py
│   ├── latex_report.py      # sezione 12, con i template TikZ della geometria
│   └── report.py            # report.md
├── tests/
└── results/
    ├── figures/
    ├── tables/
    ├── tf/                  # .npz, .mat, .json
    ├── latex/
    └── report.md
```

**Requisiti tecnici.**

- Python 3.10 o successivo, type hints ovunque, dataclass immutabili per i parametri.
- `python-control` ammesso solo come dipendenza facoltativa.
- Gli indici dei nodi sono in base 1 in tutti gli input, output, figure e documenti.
- `report.md` contiene: parametri usati con i segnaposto evidenziati, scomposizione per elemento, validazione della riduzione, frequenze e zeri, elenco delle figure, stato delle domande aperte della sezione 3 e la nota di realismo della sezione 0.
